/*
 * potential.h
 *
 *  Created on: Apr 12, 2023
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/// \file potential.h
/// \brief A table based potential for credal networks
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_POTENTIAL_H_
#define IBM_LOOPY_POTENTIAL_H_


#include <float.h>

#include "factor.h"

namespace merlin {

///
/// \brief Potential for credal networks.
///
/// Table based representation of a potential for credal networks. A potential 
/// is a set of factors defined on the same scope (set of variables) and
/// representing the convex hull (or extension) of a credal set (e.g., interval).
/// The corresponding factor table is indexed as follows (internally):
///
/// 0: [0,0,0]    4: [0,0,1]
/// 1: [1,0,0]    5: [1,0,1]
/// 2: [0,1,0]    6: [0,1,1]
/// 3: [1,1,0]    7: [1,1,1]
///
/// (For converting between source and target orders, use the convert_index class)
///
class potential {
public:
	// Typedefs:

	typedef double value;	                    ///< A factor value.
	typedef variable_set::vindex vindex;		///< Variable identifiers (0...N-1)
	typedef variable_set::vsize vsize;    		///< Variable values (0...K-1)

	// Constructors and destructor:

    ///
    /// \brief Default constructor.
    ///
    potential() {
        // empty potential
    }
    
    /// @brief Creates a scalar potential
    /// @param s a real-valued scalar (default is 1.0)
    potential(value s) {
        p_.push_back(factor(s));
        q_.push_back(factor(s));
    };

	///
	/// \brief Copy-constructor.
	///
	/// Constructs a copy from an object of the same type.
	///
	potential(potential const& f) :
			v_(f.v_), p_(f.p_), q_(f.q_) {
	};

	///
	/// \brief Class destructor
	///
	virtual ~potential() {};

	// Assignments & copy constructors:

	///
	/// \brief Assignment operator.
	///
	/// The operator performs a deep copy of the object.
	///	\param rhs	A factor object of the same type
	///	\return a reference to the current factor whose content was copied from 
	/// the factor received as argument.
	///
	potential& operator=(potential const& rhs) {
		if (this != &rhs) {
			v_ = rhs.v_;
			p_ = rhs.p_;
			q_ = rhs.q_;
		}
		return *this;
	};

	///
	/// \brief Set the domain sizes of the scope's variables.
 	///
	void set_dims() {
	};

	// Accessor functions:

	///
	/// \brief Size of the factor's scope.
	///
	/// \return the number of variables in the factor's scope.
	///
	size_t nvar() const {
		return v_.nvar();
	};
	
	///
	/// \brief Scope of the factor.
	///
	/// \return the scope (as set of variable ids) of the factor.
	///
	const variable_set& vars() const {
		return v_;
	};

	///
	/// \brief Scope of the factor.
	///
	/// \return the scope (as set of variable ids) of the factor.
	///
	const variable_set& variables() const {
		return v_;
	};
	
	///
	/// \brief Variable dimensions.
	///
	/// \return the domain sizes (i.e., dimensions) of the variables in 
	/// the factor's scope.
	///
	const vsize* dims() const {
		return v_.dims();
	};

	///
	/// \brief List of factors (the p-component).
	///
	/// \return the list of factors.
	///
	const std::vector<factor>& p() const {
		return p_;
	};


	///
	/// \brief List of factors (the q-component).
	///
	/// \return the list of factors.
	///
	const std::vector<factor>& q() const {
		return q_;
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t num_states() const {
		if (p_.empty()) {
            return 0;
        } else {
            return p_[0].num_states();
        }
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t numel() const {
        if (p_.empty()) {
            return 0;
        } else {
		    return p_[0].numel();
        }
	};

    ///
    /// \brief Add new factor to p-component
    ///
    void add_p(const factor& f) {
        p_.push_back(f); 
        v_ |= f.vars(); // update the scope
    }

    ///
    /// \brief Add new factor to q-component
    ///
    void add_q(const factor& f) {
        q_.push_back(f); 
        v_ |= f.vars(); // update the scope
    }

    ///
    /// \brief Remove dominated factors by max (p-component)
    ///
    void maximize() {
        std::vector<factor> cleaned;
        for (size_t i = 0; i < p_.size(); ++i) {
            bool found_maximal = false;
            for (size_t j = 0; j < p_.size(); ++j) {
                if (i != j && p_[j] > p_[i]) {
                    found_maximal = true;
                    break;
                }
            }

            if (!found_maximal) {
                bool found = false;
                for (size_t j = 0; j < cleaned.size(); ++j) {
                    if (p_[i] == cleaned[j]) {
                        found = true;
                    }
                }

                if (!found) {
                    cleaned.push_back(p_[i]);
                }
            }
        }

        // replace the potential's factors with the maximal ones
        p_ = cleaned;
    }

    ///
    /// \brief Remove dominated factors by max (p-component, q-component)
    ///
    void maximize2() {
		assert(p_.size() == q_.size());
        size_t num_elems = p_.size();
		std::vector<std::pair<factor, factor> > cleaned, phi;
		for (size_t i = 0; i < num_elems; ++i) {
			phi.push_back(std::make_pair(p_[i], q_[i]));
		}

        for (size_t i = 0; i < num_elems; ++i) {
            bool found_maximal = false;
            for (size_t j = 0; j < num_elems; ++j) {
                if (i != j && phi[j].first > phi[i].first && phi[j].second < phi[i].second) {
                    found_maximal = true;
                    break;
                }
            }

            if (!found_maximal) {
                bool found = false;
                for (size_t j = 0; j < cleaned.size(); ++j) {
                    if (phi[i].first == cleaned[j].first && phi[i].second == cleaned[j].second) {
                        found = true;
                    }
                }

                if (!found) {
                    cleaned.push_back(phi[i]);
                }
            }
        }

        // replace the potential's factors with the maximal ones
		p_.clear();
		q_.clear();
		for (size_t i = 0; i < cleaned.size(); ++i) {
			p_.push_back(cleaned[i].first);
			q_.push_back(cleaned[i].second);
		}
    }

    ///
    /// \brief Remove dominated factors by min (p-component)
    ///
    void minimize() {		
        std::vector<factor> cleaned;
        for (size_t i = 0; i < p_.size(); ++i) {
            bool found_minimal = false;
            for (size_t j = 0; j < p_.size(); ++j) {
                if (i != j && p_[j] < p_[i]) {
                    found_minimal = true;
                    break;
                }
            }

            if (!found_minimal) {
                bool found = false;
                for (size_t j = 0; j < cleaned.size(); ++j) {
                    if (p_[i] == cleaned[j]) {
                        found = true;
                    }
                }

                if (!found) {
                    cleaned.push_back(p_[i]);
                }
            }
        }

        // replace the potential's factors with the minimal ones
        p_ = cleaned;
    }

    ///
    /// \brief Remove dominated factors by min (p-component, q-component)
    ///
    void minimize2() {
		assert(p_.size() == q_.size());
        size_t num_elems = p_.size();
		std::vector<std::pair<factor, factor> > cleaned, phi;
		for (size_t i = 0; i < num_elems; ++i) {
			phi.push_back(std::make_pair(p_[i], q_[i]));
		}

       for (size_t i = 0; i < num_elems; ++i) {
            bool found_minimal = false;
            for (size_t j = 0; j < num_elems; ++j) {
                if (i != j && phi[j].first < phi[i].first && phi[j].second > phi[i].second) {
                    found_minimal = true;
                    break;
                }
            }

            if (!found_minimal) {
                bool found = false;
                for (size_t j = 0; j < cleaned.size(); ++j) {
                    if (phi[i].first == cleaned[j].first && phi[i].second == cleaned[j].second) {
                        found = true;
                    }
                }

                if (!found) {
                    cleaned.push_back(phi[i]);
                }
            }
        }

        // replace the potential's factors with the maximal ones
		p_.clear();
		q_.clear();
		for (size_t i = 0; i < cleaned.size(); ++i) {
			p_.push_back(cleaned[i].first);
			q_.push_back(cleaned[i].second);
		}
     }

    // Combination and marginalization operations (in place)

    ///
    /// \brief Eliminate a variable by maximization
    ///
    void max(variable v) {
        v_ /= v;
        std::vector<factor> temp;
        for (size_t i = 0; i < p_.size(); ++i) {
            factor t = p_[i].max(variable_set(v));
            temp.push_back(t);
        }

        // Replace with the new marginalized factors
        p_ = temp;
    }

    ///
    /// \brief Eliminate a variable by summation
    ///
    void sum(variable v) {
        v_ /= v;
        std::vector<factor> temp;
        for (size_t i = 0; i < p_.size(); ++i) {
            factor t = p_[i].sum(variable_set(v));
            temp.push_back(t);
        }

        // Replace with the new marginalized factors
        p_ = temp;
    }

    ///
    /// \brief Eliminate a variable by summation
    ///
    void sum2(variable v) {
        v_ /= v;
        std::vector<factor> temp1, temp2;
        for (size_t i = 0; i < p_.size(); ++i) {
            factor t = p_[i].sum(variable_set(v));
            temp1.push_back(t);
        }

        for (size_t i = 0; i < q_.size(); ++i) {
            factor t = q_[i].sum(variable_set(v));
            temp2.push_back(t);
        }

        // Replace with the new marginalized factors
        p_ = temp1;
		q_ = temp2;
    }

    ///
    /// \brief Combine two potentials by multiplication
    ///
    void multiply(const potential& B) {
        v_ |= B.vars(); // update the potential scope
        std::vector<factor> temp;
        for (size_t i = 0; i < p_.size(); ++i) {
            for (size_t j = 0; j < B.p_.size(); ++j) {
                temp.push_back(p_[i]*B.p_[j]);
            }
        }

        // Replace with the new combined factors
        p_ = temp;
    }

    ///
    /// \brief Combine two potentials by multiplication
    ///
    void multiply2(const potential& B) {
        v_ |= B.vars(); // update the potential scope
        std::vector<factor> temp1, temp2;
        for (size_t i = 0; i < p_.size(); ++i) {
            for (size_t j = 0; j < B.p_.size(); ++j) {
                temp1.push_back(p_[i]*B.p_[j]);
            }
        }
		for (size_t i = 0; i < q_.size(); ++i) {
			for (size_t j = 0; j < B.q_.size(); ++j) {
				temp2.push_back(q_[i]*B.q_[j]);
			}
		}

        // Replace with the new combined factors
        p_ = temp1;
		q_ = temp2;
    }

    /// @brief Project the potential on a variable configuration (in place)
    /// @param config current variable assignment
    void substitute(const std::map<size_t, size_t>& config) {
        std::vector<factor> temp;
        v_.clear();
		for (size_t i = 0; i < p_.size(); ++i) {
            factor f = p_[i].substitute(config);
            temp.push_back(f);
			if (v_.size() == 0) {
				v_ = f.vars(); // update the potential's scope
			}
        }

        // Replace with the new conditioned factors
        p_ = temp;
    }

    /// @brief Compute the argmax of a single variable potential
    /// @return the value that maximizes the potential
    size_t argmax() {
		// Find the most frequent value that maximizes the potential
		factor f(v_, 0.0);
		for (size_t i = 0; i < p_.size(); ++i) {
			size_t val = p_[i].argmax(); // assume factor's scope has one variable
			f[val] = f[val] + 1;
		}
		std::cout << "[ARGMAX] " << f << std::endl;
		return f.argmax();
    }

	// Boolean checks on potential properties:

	///
	/// \brief Empty potential.
	///
	/// Check if the table of factor values is empty.
	/// \return *true* if the table is empty and *false* otherwise.
	///
	bool isempty() const {
		return p_.empty();
	};

    ///
    /// @brief Check is the potential is a scalar one.
    /// @return true if the potential is a constant and false otherwise
    ///
    bool isscalar() const {
        if (v_.size() == 0) {
            if (p_.size() == 1 && p_[0].isscalar()) {
                return true;
            }
        }

        return false;
    }

	// Direct value accessor:

	///
	/// \brief Access table elements (const)
	///
	/// Direct access (read-only) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the table element at the given index.
	///		
	const factor& operator[](vsize v) const {
		return p_[v];
	};

	///
	/// \brief Rewrite table elements (non-const)
	///
	/// Direct access (read and write) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the non-const reference to table element at the given index.
	///		
	factor& operator[](vsize v) {
		return p_[v];
	};

	///
	/// \brief Access table elements (safe)
	///
	/// Direct access (read-only) to a factor value.
	/// \param i 	Index of the table element.
	/// \return the table element at the given index.
	///			
	const factor& get_p(vsize i) const {
		return p_.at(i);
	};
	const factor& get_q(vsize i) const {
		return q_.at(i);
	};


	/// @brief Get the size of the potential
	/// @return the number of elements in the potential
	size_t size() const {
		return p_.size();
	}

	///
	/// \brief Rewrite table elements (safe)
	///
	/// Re-write a factor value.
	/// \param i 	Index of the table element.
	/// \param v 	New value to be written in the table.
	///		
	void set_p(vsize i, const factor& v) {
		p_.at(i) = v;
	}
	void set_q(vsize i, const factor& v) {
		q_.at(i) = v;
	}


	///
	/// \brief Output operator (friend).
	///
	/// Write the (formatted) content of the factor to an output stream.
	/// \param out 		The output stream
	/// \param f 		The Factor to be written out
	/// \return a reference to the modified output stream containing the 
	///		content of the factor received as input.
	///
	friend std::ostream& operator<<(std::ostream& out, const potential& f) {
		out << "Potential over " << f.variables() << " is:" << std::endl;
		out << "- P components:" << std::endl;
		for (size_t j = 0; j < f.p_.size(); j++) {
			out << "  p" << j << ": " << f.p_[j] << std::endl;
		}
		out << "- Q components:" << std::endl;
		for (size_t j = 0; j < f.q_.size(); j++) {
			out << "  q" << j << ": " << f.q_[j] << std::endl;
		}

		return out;
	};


protected:

	variable_set v_;					///< Variable list vector (*scope*).
	std::vector<factor> p_;				///< List of factors (the p-component)
    std::vector<factor> q_;             ///< List of factors (the q-component)

};


inline bool compare_potentials_asc(const potential& a, const potential& b) {
	return (a.nvar() < b.nvar());
}

inline bool compare_potentials_desc(const potential& a, const potential& b) {
	return (a.nvar() > b.nvar());
}


} // namespace

#endif /* IBM_MMAP2U_POTENTIAL_H_ */
