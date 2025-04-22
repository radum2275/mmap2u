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

#ifndef IBM_MERLIN_MAP2U_H_
#define IBM_MERLIN_MAP2U_H_

#include "credal_net.h"
#include "algorithm.h"
#include "loopy2u.h"
#include "potential.h"
#include "bucket.h"
#include "search_node.h"
#include "search_space.h"
#include "bound_propagator.h"
#include "pseudotree.h"

namespace merlin {

#define SEARCH_TIMEOUT 1000


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
	
	typedef std::unique_ptr<search_node> search_node_ptr;
	typedef std::unique_ptr<bound_propagator> bound_propagator_ptr;

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
	MER_ENUM( Property , SearchMethod,PotentialApprox,Epsilon,PotentialSize,Verbose,DoCaching,DoPruning,Seed,QueryType,TimeLimit,IBound,Iterations,DoMatch,DoAndOr );


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
			set_properties("SearchMethod=bnb,PotentialApprox=none,PotentialSize=0,Epsilon=0.1,Verbose=1,DoCaching=0,DoPruning=1,Seed=0,QueryType=maximin,TimeLimit=-1,IBound=2,Iterations=1,DoMatch=0,DoAndOr=0");
			return;
		}
		m_verbose = 1;
		m_solved = false;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::SearchMethod:
				m_search_method = asgn[1]; // dfs, bb, aobb
				break;
			case Property::Epsilon:
				m_epsilon = atof(asgn[1].c_str());
				break;
			case Property::PotentialApprox:
				if (asgn[1].compare("none") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_NONE;
				} else if (asgn[1].compare("covering") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_COVERING;
				} else if (asgn[1].compare("covbound") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_COVERING_BOUND;
				} else if (asgn[1].compare("plub") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_LEAST_UPBO;
				} else if (asgn[1].compare("pglb") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_GREATEST_LOBO;
				} else if (asgn[1].compare("kmeans") == 0) {
					m_potential_approx = MERLIN_POTENTIAL_APPROX_KMEANS_BOUND;
				} else {
					std::cout << "Unsupported potential approximation scheme!" << std::endl;
				}
				break;
			case Property::PotentialSize:
				m_potential_size = atol(asgn[1].c_str());
				break;
			case Property::Iterations:
				m_iterations = atol(asgn[1].c_str());
				break;
			case Property::DoMatch:
				m_matching_strategy = atol(asgn[1].c_str());
				break;
			case Property::DoAndOr:
				m_ao_search = atol(asgn[1].c_str());
				break;
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
				break;
			case Property::DoCaching:
				m_caching = atol(asgn[1].c_str());
				break;
			case Property::DoPruning:
				m_pruning = atol(asgn[1].c_str());
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
	/// \brief Build the weighted mini-bucket heuristics
	///
	double build_heuristic();

	/// @brief Get the heuristic value corresponding to a variable
	/// @param var is the index of the variable
	/// @param assignment is the current variable assignment
	/// @return a real value representing the heuristic value
	double get_heuristic(size_t var, std::map<size_t, size_t>& assignment, bool upper);

	///
	/// \brief Moment-matching (max) in a mini-buckets partition
	///
	void moment_matching(variable vx, std::vector<potential>& partition);
	
	///
	/// \brief Max marginals
	///
	// factor maxmarginal(const factor& f, const variable_set& vs) {
	// 	return f.maxmarginal(vs);
	// }

	potential maxmarginal(const potential& p, const variable_set& vs) {
		potential result;
		for (size_t i = 0; i < p.size(); ++i) {
			factor f = p[i].maxmarginal(vs);
			result.add_p(f);
		}

		return result;
	}

	search_node* next_leaf();
	search_node* next_node();
	bool do_process(search_node* n);
	bool do_caching(search_node* n);
	bool do_pruning(search_node* n);
	bool do_expand(search_node* n);
	bool can_prune(search_node* n);
	bool generate_children(search_node* n, std::vector<search_node*>& chi);
	double heuristic(search_node* n);
	void set_cache_context(search_node* n, const std::set<size_t>& ctxt) const;
	search_node* init_search_space(double global_bound, double global_constant = 1.0);

protected:
	// Members:

	std::vector<size_t> m_order;					///< Variable elimination order
	std::map<size_t, size_t> m_evidence;			///< Evidence
	std::vector<size_t> m_query;					///< Query
	std::string m_search_method;					///< Search method (dfs, bnb, aobb, bfs, aobf)
	size_t m_verbose;								///< Verbosity level
	size_t m_seed;									///< Random number generator seed
	std::vector<int> m_best_config;					///< Best MAP config
	double m_best_cost;								///< Cost of the best MAP config
	double m_threshold;								///< Threshold used for numerical precision
	size_t m_query_type;							///< MAP type (maximin, maximax)
	double m_time_limit;							///< Time limit (default -1)
	int m_ibound;									///< Mini-buckets ibound
	double m_epsilon;								///< Epsilon value for e-covering
	size_t m_potential_size;						///< Max potential size (0 - no bounds)
	size_t m_potential_approx;						///< Potential approximation method (none, covering, lub, glb)
	size_t m_iterations;							///< Number of iterations for moment-matching
	size_t m_matching_strategy;						///< Do moment matching (0 - none, 1 - single PLUB/PGLB, 2 - exhaustive, 3 - ...)
	bool m_caching;									///< Do caching
	bool m_pruning;									///< Do pruning
	bool m_ao_search;								///< Do AND/OR search

	std::stack<search_node*> m_stack; 				///< Search stack
	std::unique_ptr<bound_propagator> m_propagator;	///< Bound propagator
	bool m_solved; 									///< Solved optimally
	size_t m_cache_hits;							///< Number of cache hits
	size_t m_num_deadends;							///< Number of deadends
	std::map<size_t, size_t> m_assignment;			///< Assignment during search
	std::unique_ptr<pseudotree> m_pseudotree;		///< Pseudo tree
	std::vector<size_t> m_domains;					///< Variable domains (including dummy)
	std::unique_ptr<search_space> m_search_space;	///< The search space

	std::vector<bucket> m_buckets;						///< The bucket structure
	std::vector<std::vector<potential>> m_intermediate;	///< The intermediate potentials
	std::vector<std::vector<potential>> m_augmented;	///< The augmented buckets
};

} // namespace




#endif /* IBM_MERLIN_MAP2U_H_ */
