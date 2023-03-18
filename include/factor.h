/*
 * factor.h
 *
 *  Created on: Feb 8, 2013
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


/// \file factor.h
/// \brief A table based factor for graphical models
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_FACTOR_H_
#define IBM_MERLIN_FACTOR_H_


#include <float.h>

#include "enum.h"
#include "utils.h"
#include "variable_set.h"
#include "index.h"

namespace merlin {

///
/// \brief Factor for graphical models.
///
/// Table based representation of a factor for graphical models. A 
/// factor encodes a potential (sometimes a probability distribution)
/// defined over a subset of discrete random variables, called a *scope*, and 
/// associates each configuration of the variables in the scope with a 
/// positive real value (sometimes a probability value). The scope is assumed
/// to be sorted lexicogaphically (e.g., [x1,x2,x3]) Also, the indexing of
/// configurations in the factor table is assumed to be based on the BigEndian
/// convention, namely the *first* variable in the ordered scope changes
/// the fastest, then the *second* variable changes its value and so on.
/// For example, consider a factor over binary variables [x1,x2,x3].
/// The corresponding factor table is indexed as follows (internally):
///
/// 0: [0,0,0]    4: [0,0,1]
/// 1: [1,0,0]    5: [1,0,1]
/// 2: [0,1,0]    6: [0,1,1]
/// 3: [1,1,0]    7: [1,1,1]
///
/// (For converting between source and target orders, use the convert_index class)
///
class factor {
public:
	// Typedefs:

	typedef std::pair<double, double> value;	///< A factor value.
	typedef variable_set::vindex vindex;		///< Variable identifiers (0...N-1)
	typedef variable_set::vsize vsize;    		///< Variable values (0...K-1)

	// Constructors and destructor:
	///
	/// \brief Copy-constructor.
	///
	/// Constructs a copy from an object of the same type.
	///
	factor(factor const& f) :
			v_(f.v_), t_(f.t_), c_(f.c_) {
	};

	///
	/// \brief Scalar constructor.
	///
	/// Creates a constant factor over an empty set of variables.
	///	\param s The scalar value used to initialize the table (default 1.0)
	///
	factor(value s = value(1.0, 1.0)) :
			v_(), t_(1, s), c_(-1) {
	};

	///
	/// \brief Constructor.
	///
	/// Creates a constant factor over a given set of variables.
	/// \param vs 	The set of variables defining the scope
	/// \param s 	The scalar value used to initialize the table (default 1.0)
	///
	factor(variable_set const& vs, value s = value(1.0, 1.0)) : v_(vs), t_(), c_(-1) {
		t_.resize(vs.num_states());
		set_dims();
		fill(s);
	};

	///
	/// \brief Constructor.
	///
	/// Creates a factor over a given set of variables and table.
	///	\param vs 	The input set of variables
	/// \param t 	The input table
	///
	factor(variable_set const& vs, value* t) :
			v_(vs), t_(), c_(-1) {
		t_.resize(v_.num_states());
		set_dims();
		std::copy(t, t + t_.size(), t_.begin());
	};

	///
	/// \brief Class destructor
	///
	virtual ~factor() {};

	// Assignments & copy constructors:

	///
	/// \brief Assignment operator.
	///
	/// The operator performs a deep copy of the object.
	///	\param rhs	A factor object of the same type
	///	\return a reference to the current factor whose content was copied from 
	/// the factor received as argument.
	///
	factor& operator=(factor const& rhs) {
		if (this != &rhs) {
			std::vector<value> tmp;
			t_.swap(tmp);                    // force vector to release memory
			v_ = rhs.v_;
			t_ = rhs.t_;
			c_ = rhs.c_;
			set_dims();                 // then reassign
		}
		return *this;
	};

	///
	/// \brief Swap the object contents.
	///
	/// Exchange the contents of "this" and the factor received as input.
	///	\param f 	The factor to exchange the content of
	///
	void swap(factor& f) {
		if (&f != this) {
			v_.swap(f.v_);
			t_.swap(f.t_);
			int tmp = c_;
			c_ = f.c_;
			f.c_ = tmp;
		}
	};

	///
	/// \brief Clone from pointer.
	///
	/// \return a pointer to the cloned factor.
	///
	virtual factor* clone() {
		factor *f = new factor(*this);
		return f;
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
	/// \brief Table of factor values.
	///
	/// \return the table containing the factor values. Each value in the table
	/// corresponds to a particular configuration of the scope variables.
	///
	const value* table() const {
		return &t_[0];
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t num_states() const {
		return t_.size();
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t numel() const {
		return t_.size();
	};

	// Boolean checks on factor properties:

	///
	/// \brief Empty factor.
	///
	/// Check if the table of factor values is empty.
	/// \return *true* if the table is empty and *false* otherwise.
	///
	bool isempty() const {
		return t_.empty();
	};

	///
	/// \brief Scalar (constant) factor.
	///
	/// Check if the factor is a *constant* (or *scalar*).
	/// \return *true* if the table size is 1 and *false* otherwise.
	///	
	bool isscalar() const {
		return numel() == 1;
	};

	// Direct value accessor:

	///
	/// \brief Access table elements (const)
	///
	/// Direct access (read-only) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the table element at the given index.
	///		
	value operator[](vsize v) const {
		return t_[v];
	};

	///
	/// \brief Rewrite table elements (non-const)
	///
	/// Direct access (read and write) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the non-const reference to table element at the given index.
	///		
	value& operator[](vsize v) {
		return t_[v];
	};

	///
	/// \brief Access table elements (safe)
	///
	/// Direct access (read-only) to a factor value.
	/// \param i 	Index of the table element.
	/// \return the table element at the given index.
	///			
	value get(vsize i) const {
		return t_.at(i);
	};

	///
	/// \brief Access the table using a scope configuration
	///
	value get_value(std::map<size_t, size_t>& config);

	///
	/// \brief Access the table using a scope configuration
	///
	void set_value(std::map<size_t, size_t>& config, value val);

	///
	/// \brief Rewrite table elements (safe)
	///
	/// Re-write a factor value.
	/// \param i 	Index of the table element.
	/// \param v 	New value to be written in the table.
	///		
	void set(vsize i, value v) {
		t_.at(i) = v;
	}

	// Filling the factor table:

	///
	/// \brief Fill with constant (in-place).
	///
	/// Set all of the factor table entries to a constant value.
	/// \param v 	The constant value to be used for filling.
	/// \return the reference to the updated factor.
	///
	factor& fill(value v) {
		std::fill(t_.begin(), t_.end(), v);
		return *this;
	};

	///
	/// \brief Fill uniformly at random (factor)
	///
	/// Create a random factor.
	/// \return the reference to the updated factor.
	///
	factor& fill_random() {
		for (size_t i = 0; i < t_.size(); ++i) {
			double v = randu();
			t_[i] = value(v, v);
		}

		return *this;
	}

	///
	/// \brief Fill a Bayesian network CPT with 1/K values.
	///
	/// Set all the factor table enties to 1/K, where K is the domain size of
	/// the child node.
	/// \return the reference to the updated factor.
	///
	factor& fill_uniform_bayes() {
		assert(c_ >= 0);
		variable_set::const_iterator i = v_.begin();
		for (; i != v_.end(); ++i) {
			if ((*i).label() == (size_t)c_) {
				double v = 1.0/(*i).states();
				std::fill(t_.begin(), t_.end(), value(v, v));
			}
		}

		return *this;
	}

	///
	/// \brief Fill a Bayesian network CPT with random values in [0, 1).
	///
	/// Set all the factor table entries to random values drawn uniformly
	/// at random between 0 and 1. The CPT is also normalized.
	/// \return the reference to the updated factor.
	///
	factor& fill_random_bayes() {
		assert(c_ >= 0);
		for (size_t i = 0; i < t_.size(); ++i) {
			double v = randu();
			t_[i] = value(v, v);
		}

		variable_set pa;
		variable ch;
		bool has_child = false;
		for (variable_set::const_iterator i = v_.begin(); i != v_.end(); ++i) {
			if (i->label() != (size_t)c_) {
				pa |= *i;
			} else {
				ch = *i;
				has_child = true;
			}
		}

		// Safety check
		assert(has_child);

		// Normalize the CPT
		if (pa.size() == 0) {
			double s = 0.0;
			for (size_t i = 0; i < t_.size(); ++i) {
				s += t_[i].first;
			}

			if (s != 0.0) {
				for (size_t i = 0; i < t_.size(); ++i) {
					t_[i].first /= s;
					t_[i].second /= s;
				}
			}
		} else {
			index_config idx1(pa); // parent configurations
			config_index idx2(v_); // factor configurations
			for (size_t i = 0; i < pa.num_states(); ++i) {

				double s = 0.0;
				std::map<size_t, size_t> cfg = idx1.convert(i);
				std::vector<size_t> indeces;
				for (size_t k = 0; k < ch.states(); ++k) {
					cfg[ch.label()] = k;
					size_t j = idx2.convert(cfg);
					indeces.push_back(j);
					s += t_[j].first;
				}

				for (size_t j = 0; j < indeces.size(); ++j) {
					size_t p = indeces[j];
					t_[p].first = (s == 0.0) ? 0.0 : (t_[p].first / s);
					t_[p].second = (s == 0.0) ? 0.0 : (t_[p].second / s);
				}
			}
		}

		return *this;
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
	friend std::ostream& operator<<(std::ostream& out, const factor& f) {
		out << "Factor over " << f.variables() << ":";
		for (size_t j = 0; j < f.t_.size(); j++) {
			out << " " << "[" << f.t_[j].first << "," << f.t_[j].second << "]";
		}
		return out;
	};

	void set_child(int c) {
		c_ = c;
	}

	int get_child() const {
		return c_;
	}

protected:

	variable_set v_;					///< Variable list vector (*scope*).
	std::vector<value> t_;				///< Table of values.
	int c_;								///< Index of the child variable (for Bayes nets only), -1 for Markov models

	///
	/// \brief Calculate the factor table size.
 	///
 	/// Compute the actual size of the factor table by multiplying the
 	/// domain sizes of the variables in its scope.
 	/// \return the factor table size.
 	///
	vsize calc_numel() const {
		vsize n = 1;
		vsize const* d = dims();
		for (size_t i = 0; i < nvar(); i++)
			n *= d[i];
		return (n > 1) ? n : 1;
	}
};


} // namespace

#endif /* IBM_MERLIN_FACTOR_H_ */
