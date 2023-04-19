/*
 * program_options.h
 *
 *  Created on: 10 September 2018
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

/// \file program_options.h
/// \brief Program options definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_PROGRAM_OPTIONS_H_
#define IBM_MERLIN_PROGRAM_OPTIONS_H_

#include <boost/program_options.hpp>

#include <iostream>
#include <sstream>

#include "base.h"

namespace po = boost::program_options;

struct ProgramOptions {

	std::string executableName;		///< Program name
	double timeLimit;				///< Time limit (in seconds)
	double memoryLimit;				///< Memory limit (in Gigs)
	int ibound;						///< Mini-bucket i-bound
	int algorithm;					///< Algorithm
	std::string scorer;				///< Scorer used for MMAP (l2u or cve2u)
	int task;						///< Inference task
	std::string modelFile;			///< Model file
	std::string evidenceFile;		///< Evidence file
	std::string queryFile;			///< Query file (MAP variables)
	std::string outputFile;			///< Output file
	size_t seed; 					///< Random number generator seed
	bool debug;						///< Debug mode
	int verbose;					///< Verbosity level (1=low, 2=medium, 3=high)
	int iterations;					///< Iterations for WMB/IJGP/JGLP
	int outputFormat;				///< Output format
	double threshold;				///< Convergence threshold
	double epsilon;					///< Epsilon value used for converting a Bayes net to an interval credal net
	double flip_probability;		///< Random flip probability (MMAP)
	std::string init_method;		///< Initialization method (MMAP)
	double init_temp;				///< Initial temperature (SA)
	double alpha;					///< Cooling factor (SA)
	size_t max_flips;				///< Max number of flips per iteration
	size_t taboo_size;				///< Max size of the taboo list 
	size_t cache_size;				///< Max size of the cache table
	size_t num_nodes;				///< Number of nodes (random problems)
	size_t num_parents;				///< Number of parents (random problems)
	size_t num_instances;			///< Number of instances (random problems)
	size_t ksize;					///< Clique size for k-trees (random problems)
	size_t kpercent;				///< Percentage of edges to be removed from k-tree (random problems)
	std::string graph_type;			///< Graph type (random problems)
	std::string query_type;			///< MMAP query type (maximin, maximax)
	size_t num_query;				///< Number of query (MAP) variables
	size_t num_samples;				///< Number of sample queries to generate
	size_t num_extras;				///< Number of extra LCN statements to generate
	size_t num_evid;				///< Number of evidence variables in LCN

public:

	// default constructor
	ProgramOptions();
};

ProgramOptions* parseCommandLine(int argc, char** argv);

inline ProgramOptions::ProgramOptions() :
		timeLimit(MERLIN_UNKNOWN),
		memoryLimit(80),
		ibound(2),
		algorithm(MERLIN_UNKNOWN),
		scorer("l2u"),
		task(MERLIN_UNKNOWN),
		seed(12345678),
		debug(false),
		verbose(0),
		iterations(10),
		outputFormat(MERLIN_OUTPUT_UAI),
		threshold(1e-06),
		epsilon(0.0),
		flip_probability(0.2),
		init_method("rand"),
		init_temp(100.0),
		alpha(0.2),
		max_flips(100),
		taboo_size(1000),
		cache_size(1000000),
		num_nodes(100),
		num_parents(3),
		num_instances(10),
		ksize(10),
		kpercent(30),
		graph_type("random"),
		query_type("maximin"),
		num_query(10),
		num_samples(10),
		num_extras(1),
		num_evid(0) {};

#endif /* IBM_MERLIN_PROGRAM_OPTIONS_H_ */
