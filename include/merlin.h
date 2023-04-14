/*
 * merlin.h
 *
 *  Created on: 20 May 2015
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

/// \file merlin.h
/// \brief Merlin inference engine interface
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef __IBM_MERLIN_H_
#define __IBM_MERLIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "credal_net.h"

///
/// Merlin probabilistic inference engine.
///
class Merlin {
	typedef size_t vindex;

protected:
	// Members:

	size_t m_task;						///< Inference task (PR, MAR, MAP, MMAP).
	size_t m_algorithm;					///< Inference algorithm
	size_t m_scorer;					///< Scorer algorithm for MMAP
	size_t m_ibound;					///< Parameter: i-bound
	size_t m_iterations;				///< Parameter: number of iterations
	size_t m_samples;					///< Parameter: number of samples
	std::string m_modelFile;			///< Model file name
	std::string m_evidenceFile;			///< Evidence file name
	std::string m_queryFile;			///< Query file name
	std::string m_outputFile;			///< Output file name
	bool m_debug;						///< Debug mode
	int m_verbose;						///< Verbosity level (default: 0)
	int m_outputFormat;					///< Output format (UAI, JSON)
	double m_threshold;					///< Threshold for early convergence
	double m_time_limit;				///< Time limit in seconds (default -1)
	double m_epsilon;					///< Epsilon value for conversion
	size_t m_seed;						///< Random number generator seed
	double m_flip_probability;			///< Random flip probability (MMAP)
	std::string m_init_method;			///< Initialization method (MMAP)
	double m_init_temp;					///< Initial temperature (SA)
	double m_alpha;						///< Cooling factor (SA)
	size_t m_max_flips;					///< Max number of flips
	size_t m_taboo_size;				///< Max size of the taboo list
	size_t m_cache_size;				///< Max size of the cache table
	size_t m_num_nodes;					///< Number of nodes (random problems)
	size_t m_num_parents;				///< Number of parents (random problems)
	size_t m_num_instances;				///< Number of instances (random problems)
	size_t m_ksize;						///< Clique size for k-trees (random problems)
	size_t m_kpercent;					///< Percentage of edges to be removed from k-tree (random problems)
	std::string m_graph_type;			///< Graph type (random problems)
	std::string m_query_type;			///< MMAP query type (maximin, maximax)
	size_t m_num_query;					///< Number of MAP query variables
	size_t m_num_samples;				///< Number of sample queries to generate
	size_t m_num_extras;				///< Number of extra LCN statements to generate
	size_t m_num_evid;					///< Number of evidence variables in LCN

private:
	// Local members:

	merlin::credal_net m_gmo;						///< Original credal net.
	std::map<vindex, size_t> m_evidence;			///< Evidence as variable value pairs.
	std::vector<vindex> m_query;					///< Query variables for MMAP tasks.
	std::string m_filename;							///< Input model file name.
	double m_ioTime;

	///
	/// \brief Clear the existing graphical model.
	///
	void clear();

	///
	/// \brief Perform safety checks.
	////
	void check();

public:

	///
	/// \brief Constructs the default Merlin engine.
	///
	Merlin();

	///
	/// \brief Destroys the Merlin engine.
	///
	~Merlin();

	///
	/// \brief Set the inference algorithm.
	/// \param alg 	The code associated with the algorithm.
	///
	void set_algorithm(size_t alg);

	///
	/// @brief Set the scorer algorithm (for MMAP)
	/// @param s the scoring method
	///
	void set_scorer(size_t s);

	///
	/// @brief Set the time limit
	/// @param t the number of seconds
	///
	void set_time_limit(double t) {
		m_time_limit = t;
	}

	///
	/// \brief Set the inference task.
	/// \param task	The code associated with the task.
	///
	void set_task(size_t task);

	///
	/// \brief Set the i-bound.
	/// \param ibound The value of the i-bound parameter.
	///
	void set_ibound(size_t ibound);

	///
	/// \brief Set the number of iterations.
	/// \param iter	The number of iterations.
	///
	void set_iterations(size_t iter);

	///
	/// \brief Set the debug mode.
	/// \param d	The flag.
	///
	void set_debug(bool v);

	///
	/// \brief Set the verbosity level.
	/// \param v	The verbosity level.
	///
	void set_verbose(int v);

	///
	/// \brief Set the convergence threshold.
	/// \param t	The threshold.
	///
	void set_threshold(double t);

	///
	/// \brief Set the epsilon threshold.
	/// \param t	The threshold.
	///
	void set_epsilon(double e);

	///
	/// \brief Set the input file name.
	/// \param f	The file name.
	///
	void set_model_file(std::string f);

	///
	/// \brief Set the output file name.
	/// \param f	The file name.
	///
	void set_output_file(std::string f);

	///
	/// \brief Set the evidence file name.
	/// \param f	The file name.
	///
	void set_evidence_file(std::string f);

	///
	/// \brief Set the query file name.
	/// \param f	The file name.
	///
	void set_query_file(std::string f);

	///
	/// \brief Set output format.
	/// \param f	The format.
	///
	void set_output_format(int f);

	///
	/// \brief Set the random number generator seed
	///
	void set_seed(size_t s);

	///
	/// \brief Set the random flip probability
	///
	void set_flip_probability(double p);

	///
	/// \brief Set the initialization method
	///
	void set_init_method(std::string m);

	///
	/// \brief Set the initial temperature
	///
	void set_init_temp(double t);

	///
	/// \brief Set the cooling factor
	///
	void set_alpha(double a);

	///
	/// \brief Set the max number of flips
	///
	void set_max_flips(size_t m);

	///
	/// \brief Set the taboo list size
	///
	void set_taboo_size(size_t m);

	///
	/// \brief Set the cache table size
	///
	void set_cache_size(size_t m);

	void set_num_nodes(size_t n);
	void set_num_parents(size_t p);
	void set_num_instances(size_t m);
	void set_ksize(size_t k);
	void set_kpercent(size_t p);
	void set_graph_type(std::string s);
	void set_query_type(std::string s);
	void set_num_query(size_t q);
	void set_num_samples(size_t s);
	void set_num_extras(size_t e);
	void set_num_evid(size_t e);

	///
	/// \brief Initialize the solver.
	///	\return *true* if succesful and *false* otherwise.
	///
	bool init();

	///
	/// \brief Solve the inference task given current evidence.
	///	\return 0 if succesful and 1 otherwise.
	///
	int run();


protected:

	///
	/// \brief Read the graphical model from a file in the UAI format.
	/// \param file_name	The input file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_model(const char* filename);

	///
	/// \brief Read the evidence.
	/// \param file_name	The evidence file name.
	///	\return *true* if successful and *false* otherwise.
	//
	bool read_evidence(const char* filename);

	///
	/// \brief Read the query variables (MMAP task and joint marginals).
	/// \param file_name	The query file name.
	///	\return *true* if successful and *false* otherwise.
	///
	bool read_query(const char* filename);

};

#endif /* __IBM_MERLIN_H_ */
