/*
 * cve2u.h
 *
 *  Created on: 04 Apr 2023
 *      Author: radu
 *
 * Copyright (c) 2023, International Business Machines Corporation
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

/// \file cve2u.h
/// \brief CVE2U for credal nets with intervals and binary variables
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_CVE2U_H_
#define IBM_LOOPY_CVE2U_H_

#include "credal_net.h"
#include "algorithm.h"
#include "bucket.h"

namespace merlin {


/**
 * Credal Variable Elimination 2U (CVE2U)
 *
 * Tasks supported: PE (probability of evidence)
 *
 */

class cve2u: public credal_net, public algorithm {
public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	cve2u() : credal_net() {
		set_properties();
	}

	///
	/// \brief Constructor with a credal net.
	///
	cve2u(const credal_net& cn) : credal_net(cn) {
		set_properties();
	}

	///
	/// \brief Constructor with a set of factors.
 	///
	cve2u(const std::vector<interval>& fs) : credal_net(fs) {
		set_properties();
	}

	///
	/// \brief Destructor
	///
	~cve2u() {
	};

	inline const interval& belief(size_t i) const {
		return m_beliefs[i];
	}
	inline const interval& belief(variable v) const {
		return m_beliefs[v];
	}
	inline const std::vector<interval>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Write the solution to the output stream.
	/// \param out		 		The output stream
	/// \param output_format	The output format (json or uai)
	///
	void write_solution(std::ostream& out, int output_format);

	///
	/// \brief Run the CVE2U algorithm for P(e).
	///
	void run();

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Verbose,Seed );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the variable order.
	///
	inline void set_order(const variable_order_t& ord) {
		m_order = ord;
	}

	///
	/// \brief Get the variable order.
	///
	inline const variable_order_t& get_order() const {
		return m_order;
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Verbose=1,Seed=0");
			return;
		}
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
				break;
			case Property::Seed:
				m_seed = atol(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Reset the internal state of the solver
	///
	void reset();
	void reset(const std::map<size_t, size_t>& config);

	///
	/// @brief Compute the probability of evidence
	/// @param config a variable assignment
	/// @return lower and upper bounds on the probability of evidence
	///
	std::pair<double, double> eval(const std::map<size_t, size_t>& config);

	///
	/// \brief Set the evidence variables
	///
	void set_evidence(const std::map<size_t, size_t>& ev) {
		m_evidence = ev;
	}

	///
	/// \brief Initialize the CVE 2U algorithm.
	///
	void init();

	///
	/// \brief Update the beliefs for all variables
	///
	void update_beliefs();

	///
	/// \brief Update the belief of a variable
	///
	void upate_belief(variable x);
		

protected:

	double ve(bool upper = true);

protected:
	// Members:

	variable_order_t m_order;						///< Variable order
	std::vector<interval> m_beliefs; 				///< Marginals
	std::map<size_t, size_t> m_evidence;			///< Evidence
	std::vector<bucket> m_buckets;					///< Bucket structure
	std::vector<size_t> m_positions;				///< Positions of variables in ordering
	size_t m_verbose;								///< Verbosity level
	size_t m_seed;									///< Random number generator seed
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
