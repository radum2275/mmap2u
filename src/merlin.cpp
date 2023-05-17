/*
 * merlin.cpp
 *
 *  Created on: 15 May 2015
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

// Merlin library core.
#include "utils.h"
#include "merlin.h"
#include "loopy2u.h"
#include "ipe2u.h"
#include "cve2u.h"
#include "bn2cn.h"
#include "mmap2u.h"
#include "generator.h"

#include <stdlib.h>

///
/// \brief Constructs the default Merlin engine.
///
Merlin::Merlin() {
	m_task = MERLIN_TASK_MAR;
	m_algorithm = MERLIN_ALGO_L2U;
	m_scorer = "l2u";
	m_ibound = 2;
	m_iterations = 10;
	m_time_limit = -1;
	m_samples = 100;
	m_debug = false;
	m_verbose = 0;
	m_outputFormat = MERLIN_OUTPUT_UAI;
	m_ioTime = 0;
	m_threshold = 1e-06;
	m_epsilon = 0.0;
	m_seed = 12345678;
	m_flip_probability = 0.2;
	m_init_method = "rand";
	m_init_temp = 100;
	m_alpha = 0.2;
	m_max_flips = 100;
	m_taboo_size = 1000;
	m_cache_size = 1000000;
	m_num_nodes = 100;
	m_num_parents = 3;
	m_num_instances = 10;
	m_ksize = 10;
	m_kpercent = 30;
	m_graph_type = "random";
	m_query_type = "maximin";
	m_num_query = 10;
	m_num_samples = 10;
	m_num_extras = 1;
	m_num_evid = 0;
}

///
/// \brief Destroys the Merlin engine.
///
Merlin::~Merlin() {
	clear();
}

///
/// \brief Clears the internal structures.
///
void Merlin::clear() {
}

///
/// \brief Set the random number generator seed
///
void Merlin::set_seed(size_t s) {
	m_seed = s;
}

///
/// \brief Set the inference algorithm.
/// \param alg 	The code associated with the algorithm.
///
void Merlin::set_algorithm(size_t alg) {
	m_algorithm = alg;
}

///
/// @brief Set the scorer algorithm (for MMAP)
/// @param s the scoring method code
///
void Merlin::set_scorer(std::string s) {
	m_scorer = s;
}

///
/// \brief Set the inference task.
/// \param task	The code associated with the task.
///
void Merlin::set_task(size_t task) {
	m_task = task;
}

///
/// \brief Set the i-bound.
/// \param ibound The value of the i-bound parameter.
///
void Merlin::set_ibound(size_t ibound) {
	m_ibound = ibound;
}

///
/// \brief Set the number of iterations.
/// \param iter	The number of iterations.
///
void Merlin::set_iterations(size_t iter) {
	m_iterations = iter;
}

///
/// \brief Set the debug flag.
///
void Merlin::set_debug(bool v) {
	m_debug = v;
}

///
/// \brief Set the verbosity level.
///
void Merlin::set_verbose(int v) {
	m_verbose = v;
}

///
/// \brief Set the convergence threshold.
///
void Merlin::set_threshold(double t) {
	m_threshold = t;
}

///
/// \brief Set the epsilon threshold.
///
void Merlin::set_epsilon(double e) {
	m_epsilon = e;
}

///
/// \brief Set the input file name.
/// \param f	The file name.
///
void Merlin::set_model_file(std::string f) {
	m_modelFile = f;
}

///
/// \brief Set the output file name.
/// \param f	The file name.
///
void Merlin::set_output_file(std::string f) {
	m_outputFile = f;
}

///
/// \brief Set the evidence file name.
/// \param f	The file name.
///
void Merlin::set_evidence_file(std::string f) {
	m_evidenceFile = f;
}

///
/// \brief Set the query file name.
/// \param f	The file name.
///
void Merlin::set_query_file(std::string f) {
	m_queryFile = f;
}

///
/// \brief Set the output format.
///
void Merlin::set_output_format(int f) {
	m_outputFormat = f;
}

///
/// \brief Set the random flip probability.
///
void Merlin::set_flip_probability(double p) {
	m_flip_probability = p;
}

///
/// \brief Set the initialization method.
///
void Merlin::set_init_method(std::string m) {
	m_init_method = m;
}

///
/// \brief Set the initial temperature
///
void Merlin::set_init_temp(double t) {
	m_init_temp = t;
}

///
/// \brief Set the cooling factor
///
void Merlin::set_alpha(double a) {
	m_alpha = a;
}

///
/// \brief Set the max number of flips
///
void Merlin::set_max_flips(size_t m) {
	m_max_flips = m;
}

///
/// \brief Set the taboo list size
///
void Merlin::set_taboo_size(size_t m) {
	m_taboo_size = m;
}

///
/// \brief Set the cache table size
///
void Merlin::set_cache_size(size_t m) {
	m_cache_size = m;
}

void Merlin::set_num_nodes(size_t n) {
	m_num_nodes = n;
}
void Merlin::set_num_parents(size_t p) {
	m_num_parents = p;
}
void Merlin::set_num_instances(size_t m) {
	m_num_instances = m;
}
void Merlin::set_ksize(size_t k) {
	m_ksize = k;
}
void Merlin::set_kpercent(size_t p) {
	m_kpercent = p;
}
void Merlin::set_graph_type(std::string s) {
	m_graph_type = s;
}
void Merlin::set_query_type(std::string s) {
	m_query_type = s;
}
void Merlin::set_num_query(size_t q) {
	m_num_query = q;
}
void Merlin::set_num_samples(size_t s) {
	m_num_samples = s;
}
void Merlin::set_num_extras(size_t e) {
	m_num_extras = e;
}
void Merlin::set_num_evid(size_t e) {
	m_num_evid = e;
}

///
/// \brief Read the credal net.
/// \param filename	The input file name.
///
bool Merlin::read_model(const char* filename) {
	try {

		// Read the credal net
		std::cout << "[MERLIN] Reading file: " << filename << std::endl;
		m_filename = std::string(filename);
		std::ifstream is(filename);
		if (is.fail()) {
			std::string err_msg("Cannot open the input file: ");
			err_msg += std::string(filename);
			throw std::runtime_error(err_msg);
		}

		m_gmo.read(is); // throws a runtime_error in case of failure

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Read the evidence as variable-value pairs.
/// \param filename	The evidence file name.
///
bool Merlin::read_evidence(const char* filename) {
	try {

		// Open the evidence file
		std::ifstream in(filename);
		if (in.fail()) {
			std::string err_msg("Cannot open the evidence file: ");
			err_msg += filename;
			throw std::runtime_error(err_msg);
		}

		// Clear any previous evidence (m_evidence is a map)
		m_evidence.clear();

		// Read the evidence file
		int num_evid;
		in >> num_evid;
		for (int i = 0; i < num_evid; ++i) {
			vindex var;
			size_t val;
			in >> var >> val;
			m_evidence[var] = val;
		}

		// Close the evidence file
		in.close();

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}


///
/// \brief Read the query variables (MMAP task only).
/// \param filename	The query file name.
///
bool Merlin::read_query(const char* filename) {
	try {

		// Open the query file
		std::ifstream in(filename);
		if (in.fail()) {
			throw std::runtime_error("Error while opening the query file.");
		}

		// Clear any previous query
		m_query.clear();

		// Read the query file
		int num_vars;
		in >> num_vars;
		for (int i = 0; i < num_vars; ++i) {
			vindex var;
			in >> var;
			m_query.push_back(var);
		}

		// Close the evidence file
		in.close();

		// Sort the query variables in ascending order
		std::sort(m_query.begin(), m_query.end());

		return true;
	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}
}

///
/// \brief Safety checks.
///
void Merlin::check() {
	// do nothing
}

///
/// \brief Initialize the solver.
///
bool Merlin::init() {

	double timestamp = merlin::timeSystem();
	if (!read_model(m_modelFile.c_str())) {
		return false;
	}
	if (m_evidenceFile.empty() == false) {
		if (!read_evidence(m_evidenceFile.c_str())) {
			return false;
		}
	}
	if (m_queryFile.empty() == false) {
		if (!read_query(m_queryFile.c_str())) {
			return false;
		}
	}

	m_ioTime = (merlin::timeSystem() - timestamp);
	merlin::rand_seed(m_seed);
	
	return true;
}


///
/// \brief Solve the inference task given current evidence.
///
int Merlin::run() {

	try {

		// Safety checks
		check();

		// Prologue
		std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;
		std::cout << "[MERLIN] Starting up the engine ..." << std::endl;

		// Set the output file
		if (m_task == MERLIN_TASK_MAR) {
			if (m_outputFile.empty()) {
				size_t found = m_filename.find_last_of("/");
				std::string prob_name = (found != std::string::npos) ?
						m_filename.substr(found + 1) : m_filename;
				m_outputFile = "./" + prob_name;
			}

			// Set the output format
			m_outputFile += ".MAR";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}
		} else if (m_task == MERLIN_TASK_MMAP) {
			if (m_outputFile.empty()) {
				size_t found = m_filename.find_last_of("/");
				std::string prob_name = (found != std::string::npos) ?
						m_filename.substr(found + 1) : m_filename;
				m_outputFile = "./" + prob_name;
			}

			// Set the output format
			m_outputFile += ".MMAP";
			if (m_outputFormat == MERLIN_OUTPUT_JSON) {
				m_outputFile += ".json";
			}
		} else if (m_task == MERLIN_TASK_CONV) {
			assert(m_outputFile.empty() == false);
		} else if (m_task == MERLIN_TASK_GEN) {
			// do nothing
		}

		// Open the output stream
		if (m_task != MERLIN_TASK_GEN) {
			std::ofstream out(m_outputFile.c_str());
			if (out.fail()) {
				std::string err_msg("Cannot open output file: ");
				err_msg += m_outputFile;
				throw std::runtime_error(err_msg);
			}
		
			if (m_algorithm == MERLIN_ALGO_L2U) {
				merlin::loopy2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.run();
				s.write_solution(out, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_IPE2U) {
				merlin::ipe2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.run();
				s.write_solution(std::cout, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_CVE2U) {
				merlin::cve2u s(m_gmo);
				std::ostringstream oss;
				oss << "Verbose=" << m_verbose << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.run();
			} else if (m_algorithm == MERLIN_ALGO_MMAP_HILL) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=hc,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_MMAP_TABOO) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=ts,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "TabooSize=" << m_taboo_size << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_MMAP_SA) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=sa,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat);
			} else if (m_algorithm == MERLIN_ALGO_MMAP_CVE) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=cve,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat); 
			} else if (m_algorithm == MERLIN_ALGO_MMAP_CMBE) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=cmbe,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "IBound=" << m_ibound << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat); 
			} else if (m_algorithm == MERLIN_ALGO_MMAP_DFS) {
				merlin::mmap2u s(m_gmo);
				std::ostringstream oss;
				oss << "StopIter=" << m_iterations << ","
					<< "FlipProb=" << m_flip_probability << ","
					<< "InitMethod=" << m_init_method << ","
					<< "InitTemp=" << m_init_temp << ","
					<< "Alpha=" << m_alpha << ","
					<< "MaxFlips=" << m_max_flips << ","
					<< "SearchMethod=naive,"
					<< "Scorer=" << m_scorer << ","
					<< "Threshold=" << m_threshold << ","
					<< "Verbose=" << m_verbose << ","
					<< "QueryType=" << m_query_type << ","
					<< "CacheSize=" << m_cache_size << ","
					<< "TimeLimit=" << m_time_limit << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.set_evidence(m_evidence);
				s.set_query(m_query);
				s.run();
				s.write_solution(std::cout, m_outputFormat); 
			} else if (m_algorithm == MERLIN_ALGO_CONVERT) {
				merlin::bn2cn s(m_gmo);
				std::ostringstream oss;
				oss << "Epsilon=" << m_epsilon << ","
					<< "Verbose=" << m_verbose << ","
					<< "Seed=" << m_seed;
				s.set_properties(oss.str());
				s.run();
				s.write_solution(out, MERLIN_OUTPUT_UAI);
			}

			out.close();
		} else { 	
			merlin::generator s;
			std::ostringstream oss;
			oss << "Epsilon=" << m_epsilon << ","
				<< "Seed=" << m_seed << ","
				<< "Graph=" << m_graph_type << ","
				<< "Instances=" << m_num_instances << ","
				<< "Nodes=" << m_num_nodes << ","
				<< "Parents=" << m_num_parents << ","
				<< "KSize=" << m_ksize << ","
				<< "KPercent=" << m_kpercent << ","
				<< "Query=" << m_num_query << ","
				<< "Samples=" << m_num_samples << ","
				<< "Extras=" << m_num_extras << ","
				<< "Evid=" << m_num_evid;
			s.set_properties(oss.str());
			s.set_input_filename(m_filename);
			s.run();
		}

		std::cout << "[MERLIN] I/O time is " << std::fixed
				<< std::setprecision(MERLIN_PRECISION) << m_ioTime
				<< " seconds" << std::endl;

		return MERLIN_EXIT_SUCCESS;

	} catch (const std::runtime_error& e) {
		std::cerr << e.what() << std::endl;
		return MERLIN_EXIT_FAILURE;
	}
}


