/*
 * mmap2u.h
 *
 *  Created on: 29 Oct 2021
 *      Author: radu
 *
 * Copyright (c) 2021, International Business Machines Corporation
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

/// \file mmap2u.h
/// \brief Hill climbing for MMAP in credal nets with intervals and binary variables
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_MMAP2U_H_
#define IBM_LOOPY_MMAP2U_H_

#include "credal_net.h"
#include "algorithm.h"
#include "loopy2u.h"

namespace merlin {


/**
 * Local search for credal MMAP
 *
 * Tasks supported: MMAP
 *
 */

class mmap2u: public credal_net, public algorithm {
public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	mmap2u() : credal_net() {
		set_properties();
		m_scorer = NULL;
	}

	///
	/// \brief Constructor with a credal net.
	///
	mmap2u(const credal_net& cn) : credal_net(cn) {
		set_properties();
		m_scorer = NULL;
	}

	///
	/// \brief Destructor
	///
	~mmap2u() {
		if (m_scorer != NULL) {
			delete m_scorer;
		}
	};

	inline const factor& belief(size_t i) const {
		return m_beliefs[i];
	}
	inline const factor& belief(variable v) const {
		return m_beliefs[v];
	}
	inline const std::vector<factor>& beliefs() const {
		return m_beliefs;
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
	MER_ENUM( Property , StopIter,FlipProb,InitMethod,InitTemp,Alpha,MaxFlips,SearchMethod,Threshold,Verbose,Seed,TabooSize,CacheSize,QueryType );


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
			set_properties("StopIter=10,FlipProb=0.2,InitTemp=100,Alpha=0.2,MaxFlips=100,InitMethod=rand,SearchMethod=hc,Threshold=1e-06,Verbose=1,Seed=0,TabooSize=100,CacheSize=100,QueryType=maximin");
			return;
		}
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::StopIter:
				m_iterations = (size_t) atol(asgn[1].c_str());
				break;
			case Property::FlipProb:
				m_flip_probability = atof(asgn[1].c_str());
				break;
			case Property::InitMethod:
				m_init_method = asgn[1]; // rand, mpe, mle
				break;
			case Property::SearchMethod:
				m_search_method = asgn[1]; // hc, ts, sa
				break;
			case Property::Alpha:
				m_alpha = atof(asgn[1].c_str());
				break;
			case Property::InitTemp:
				m_init_temperature = atof(asgn[1].c_str());
				break;
			case Property::MaxFlips:
				m_max_flips = (size_t) atol(asgn[1].c_str());
				break;
			case Property::TabooSize:
				m_taboo_size = (size_t) atol(asgn[1].c_str());
				break;
			case Property::CacheSize:
				m_cache_size = (size_t) atol(asgn[1].c_str());
				break;
			case Property::Threshold:
				m_threshold = atof(asgn[1].c_str());
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
				break;
			case Property::Seed:
				m_seed = atol(asgn[1].c_str());
				break;
			case Property::QueryType:
				if (asgn[1].compare("maximax") == 0) {
					m_query_type = MERLIN_MMAP_MAXIMAX;
				} else if (asgn[1].compare("maximin") == 0) {
					m_query_type = MERLIN_MMAP_MAXIMIN;
				} else {
					// TODO: need to handle non-dominated intervals
					m_query_type = MERLIN_MMAP_INTERVAL;
				}
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
	/// \brief Set the query variables (e.g., MAP variables for MMAP)
	///
	void set_query(const std::vector<size_t>& query) {
		m_query = query;
	}

protected:

	///
	/// \brief Hill climbing search
	///
	void hill_climbing();

	///
	/// \brief Taboo search
	///
	void taboo_search();
	
	///
	/// \brief Simulated annealing
	///
	void simulated_annealing();
	
	///
	/// \brief Calculate the score of a MAP configuration
	///
	double score(const std::vector<size_t>& config);

	///
	/// \brief Compute initial MAP configuration
	///
	std::vector<size_t> init_config();

	///
	/// \brief Get a random neighbor of the input MAP configuration
	///
	std::vector<size_t> get_random_neighbor(const std::vector<size_t>& config);

	///
	/// \brief Find the neighborhood of the input MAP configuration
	///
	void find_neighbors(variable x, const std::vector<size_t>& config, 
		std::vector<std::vector<size_t> >& neighbors);

	///
	/// \brief Find the neighborhood of the input MAP configuration
	///
	void find_neighbors(const std::vector<size_t>& config, 
    	std::vector<std::vector<size_t> >& neighbors);

	///
	/// \brief Create a string key from a variable assignment
	///
	std::string make_key(const std::vector<size_t>& config);
	
	///
	/// \brief Convert a variable assignment to a string
	///
	std::string to_string(variable_set& vars, std::map<size_t, size_t>& config);

	///
	/// \brief Initialize scorer
	///
	void init_scorer();

	///
	/// \brief Update the auxiliary factors
	///
	void update_aux_factors(const std::map<size_t, size_t>& config);

protected:
	// Members:

	variable_order_t m_order;						///< Variable order
	std::vector<factor> m_beliefs; 					///< Marginals
	std::map<size_t, size_t> m_evidence;			///< Evidence
	std::vector<size_t> m_query;					///< Query
	size_t m_iterations;							///< Number of iterations
	double m_threshold;								///< Convergence threshold (default=1e-06)
	double m_flip_probability;						///< Random flip probability
	std::string m_init_method;						///< Initialization method (rand, mpe, mle)
	std::string m_search_method;					///< Search method (hill, taboo, aneal)
	std::vector<size_t> m_schedule;					///< Propagation schedule
	size_t m_verbose;								///< Verbosity level
	size_t m_seed;									///< Random number generator seed
	std::vector<size_t> m_best_config;				///< Best MMAP config
	double m_best_score;							///< Best MMAP score
	double m_init_temperature;						///< [SA] Initial temperature
	double m_alpha;									///< [SA] Cooling factor
	size_t m_max_flips;								///< Max number of flips
	size_t m_taboo_size;							///< Max size of the taboo list
	size_t m_cache_size;							///< Max size of the cache table
	size_t m_query_type;							///< MMAP type (maximin, maximax, interval)
	loopy2u* m_scorer;								///< Scorer
	std::vector<findex> m_aux_fid;					///< List of auxiliary factor indices
	std::vector<vindex> m_aux_vid;					///< List of auxiliary variable indices
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
