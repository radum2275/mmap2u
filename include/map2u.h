/*
 * map2u.h
 *
 *  Created on: 18 Dec 2024
 *      Author: radu
 *
 * Copyright (c) 2024, International Business Machines Corporation. All rights reserved.
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

/// \file map2u.h
/// \brief Branch and Bound algorithms for MAP in credal nets with intervals and binary variables
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_MAP2U_H_
#define IBM_LOOPY_MAP2U_H_

#include "credal_net.h"
#include "algorithm.h"
#include "loopy2u.h"
#include "potential.h"
#include "bucket.h"

namespace merlin {


/**
 * Branch and Bound search algorithms for MAP in credal networks.
 *
 * Tasks supported: MAP
 *
 */

class map2u: public credal_net, public algorithm {
public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	map2u() : credal_net() {
		set_properties();
	}

	///
	/// \brief Constructor with a credal net.
	///
	map2u(const credal_net& cn) : credal_net(cn) {
		set_properties();
	}

	///
	/// \brief Destructor
	///
	~map2u() {
	};

	inline const interval& belief(size_t i) const {
		throw 0;
	}
	inline const interval& belief(variable v) const {
		throw 0;
	}
	inline const std::vector<interval>& beliefs() const {
		throw 0;
	}

	///
	/// \brief Write the solution to the output stream.
	/// \param out		 		The output stream
	/// \param output_format	The output format (json or uai)
	///
	void write_solution(std::ostream& out, int output_format);

	///
	/// \brief Initialize the credal MMAP algorithm.
	///
	void init();

	///
	/// \brief Run the credal MMAP algorithm.
	///
	void run();

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , SearchMethod,PotentialApprox,Epsilon,PotentialSize,Verbose,Seed,QueryType,TimeLimit,IBound );


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
			set_properties("SearchMethod=bnb,PotentialApprox=none,PotentialSize=0,Epsilon=0.1,Verbose=1,Seed=0,QueryType=maximin,TimeLimit=-1,IBound=2");
			return;
		}
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::SearchMethod:
				m_search_method = asgn[1]; // hc, ts, sa
				break;
			case Property::Epsilon:
				m_epsilon = atof(asgn[1].c_str());
				break;
			case Property::PotentialApprox:
				if (asgn[1].compare("none") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_NONE;
				} else if (asgn[1].compare("covering") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_COVERING;
				} else if (asgn[1].compare("lpub") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_LEAST_UPBO;
				} else if (asgn[1].compare("gplb") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_GREATEST_LOBO;
				} else {
					std::cout << "Unsupported potential approximation scheme!" << std::endl;
				}
				break;
			case Property::PotentialSize:
				m_potential_size = atol(asgn[1].c_str());
				break;
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
				break;
			case Property::Seed:
				m_seed = atol(asgn[1].c_str());
				break;
			case Property::QueryType:
				if (asgn[1].compare("maximax") == 0) {
					m_query_type = MERLIN_MAP_MAXIMAX;
				} else if (asgn[1].compare("maximin") == 0) {
					m_query_type = MERLIN_MAP_MAXIMIN;
				} else {
					std::cout << "Only maximin and maximax MAP queries are supported!" << std::endl;
				}
				break;
			case Property::TimeLimit:
				m_time_limit = atof(asgn[1].c_str());
				break;
			case Property::IBound:
				m_ibound = atoi(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Set the evidence variables
	///
	void set_evidence(const std::map<size_t, size_t>& ev) {
		m_evidence = ev;
	}

	///
	/// \brief Set the query variables (e.g., MAP variables)
	///
	void set_query(const std::vector<size_t>& query) {
		m_query = query;
	}

protected:

	///
	/// \brief OR Branch and Bound Search
	///
	void bnb();

	///
	/// \brief Depth-First search
	///
	void dfs();
	
	///
	/// \brief AND/OR Branch and Bound Search
	///
	void aobb();
	
	///
	/// \brief Weighted Mini-Buckets
	///
	void wmb();

	///
	/// \brief Convert a variable assignment to a string
	///
	std::string to_string(variable_set& vars, std::map<size_t, size_t>& config);

	///
	/// \brief Precompile the weighted mini-bucket heuristics for MAP
	///
	void precompile_heuristics();

protected:
	// Members:

	variable_order_t m_order;						///< Variable elimination order
	std::map<size_t, size_t> m_evidence;			///< Evidence
	std::vector<size_t> m_query;					///< Query
	std::string m_search_method;					///< Search method (dfs, bnb, aobb, bfs, aobf)
	size_t m_verbose;								///< Verbosity level
	size_t m_seed;									///< Random number generator seed
	std::vector<int> m_best_config;					///< Best MAP config
	double m_best_score;							///< Score of the best MAP config
	double m_threshold;								///< Threshold used for numerical precision
	size_t m_query_type;							///< MAP type (maximin, maximax)
	double m_time_limit;							///< Time limit (default -1)
	int m_ibound;									///< Mini-buckets ibound
	double m_epsilon;								///< Epsilon value for e-covering
	size_t m_potential_size;						///< Max potential size (0 - no bounds)
	size_t m_potential_approx;						///< Potential approximation method (none, covering, lub, glb)
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
