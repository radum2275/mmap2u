/*
 * program_options.cpp
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

/// \file program_options.cpp
/// \brief Program options definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#include "utils.h"
#include "program_options.h"
#include <ostream>

///
/// \brief Parse the command line arguments.
///
ProgramOptions* parseCommandLine(int argc, char** argv) {

	ProgramOptions* opt = new ProgramOptions;

	// executable name
	opt->executableName = argv[0];

	try {
		
		po::options_description desc("Valid options");
		desc.add_options()
			("input-file,f", po::value<std::string>(), "path to problem file (required)")
			("evidence-file,e", po::value<std::string>(), "path to evidence file (required)")
			("query-file,q", po::value<std::string>(), "path to query file file (optional)")
			("output-file,o", po::value<std::string>(), "path to output file (optional)")
			("algorithm,a", po::value<std::string>(), "inference algorithm (required): bte, cte, wmb, ijgp, lbp, jglp, gibbs")
			("task,t", po::value<std::string>(), "inference task (use PR, MAR, MAP, MMAP)")
			("time-limit,l", po::value<int>(), "time limit in seconds")
			("seed,s", po::value<size_t>(), "seed for the random number generator")
			("verbose,v", po::value<int>(), "specify verbosity level")
			("debug,d", "enable debug mode")
			("iterations,n", po::value<int>(), "number of iterations")
			("threshold,T", po::value<double>(), "threshold for L2U convergence")
			("epsilon", po::value<double>(), "epsilon for converting to an interval credal net")
			("flip-proba", po::value<double>(), "random flip probability for MMAP")
			("init-method", po::value<std::string>(), "initialization method for MMAP")
			("output-format,O", po::value<std::string>(), "output file format (required)")
			("alpha", po::value<double>(), "cooling factor for SA")
			("init-temp", po::value<double>(), "initial temperature for SA")
			("max-flips", po::value<size_t>(), "max number of flips per iteration")
			("taboo-size", po::value<size_t>(), "max configurations in the taboo list")
			("cache-size", po::value<size_t>(), "max configurations in the cache")
			("nodes", po::value<size_t>(), "number of nodes")
			("parents", po::value<size_t>(), "number of parents")
			("instances", po::value<size_t>(), "number of instances")
			("ksize", po::value<size_t>(), "clique size for k-trees")
			("kpercent", po::value<size_t>(), "percentage of edges to be removed from k-tree")
			("graph-type", po::value<std::string>(), "graph type for random problems: random, grid, ktree")
			("query-type", po::value<std::string>(), "MMAP query type: maximin, maximax")			
			("num-query", po::value<size_t>(), "number of MAP query variables")
			("num-samples", po::value<size_t>(), "number of sample queries to generate")
			("num-extras", po::value<size_t>(), "number of extra LCN statements to generate")
			("num-evid", po::value<size_t>(), "number of extra LCN evidence sentences")
			("help,h", "produces this help message");

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		// parse help
		if (vm.count("help")) {
			std::cout << std::endl << desc << std::endl;
			delete opt;
			exit(0);
		}

		// parse verbosity level
		if (vm.count("verbose")) {
			opt->verbose = vm["verbose"].as<int>();
		}

		// parse convergence threshold
		if (vm.count("threshold")) {
			opt->threshold = vm["threshold"].as<double>();
		}

		// parse debug mode (switch)
		if (vm.count("debug")) {
			opt->debug = true;
		}

		// parse input file
		// if (!vm.count("input-file")) {
		// 	std::string err_msg("Input model file is required. ");
		// 	err_msg += "Call with '" + std::string(argv[0]) + " --help' ";
		// 	err_msg += "for full description of the command line arguments.";
		// 	throw std::runtime_error(err_msg);
		// } else {
		// 	opt->modelFile = vm["input-file"].as<std::string>();
		// }
		if (vm.count("input-file")) {
			opt->modelFile = vm["input-file"].as<std::string>();
		}

		// parse the query variables file
		if (vm.count("query-file")) {
			opt->queryFile = vm["query-file"].as<std::string>();
		}

		// parse the evidence file
		if (vm.count("evidence-file")) {
			opt->evidenceFile = vm["evidence-file"].as<std::string>();
		}

		// parse output file
		if (vm.count("output-file")) {
			opt->outputFile = vm["output-file"].as<std::string>();
		}

		// parse inference task
		if (vm.count("task")) {
			std::string task = vm["task"].as<std::string>();
			if (task.compare("MAR") == 0) {
				opt->task = MERLIN_TASK_MAR;
			} else if (task.compare("MMAP") == 0) {
				opt->task = MERLIN_TASK_MMAP;
			} else if (task.compare("CONV") == 0) {
				opt->task = MERLIN_TASK_CONV;
			} else if (task.compare("GEN") == 0) {
				opt->task = MERLIN_TASK_GEN;
			} else {
				std::string err_msg("Inference task ");
				err_msg += task + " is not supported.";
				throw std::runtime_error(err_msg);
			}
		}

		// parse algorithm
		if (vm.count("algorithm")) {
			std::string alg = vm["algorithm"].as<std::string>();
			if (alg.compare("l2u") == 0) {
				opt->algorithm = MERLIN_ALGO_L2U;
			} else if (alg.compare("ipe2u") == 0) {
				opt->algorithm = MERLIN_ALGO_IPE2U;
			} else if (alg.compare("sv2u") == 0) {
				opt->algorithm = MERLIN_ALGO_SV2U;
			} else if (alg.compare("hc") == 0) {
				opt->algorithm = MERLIN_ALGO_MMAP_HILL;
			} else if (alg.compare("ts") == 0) {
				opt->algorithm = MERLIN_ALGO_MMAP_TABOO;
			} else if (alg.compare("sa") == 0) {
				opt->algorithm = MERLIN_ALGO_MMAP_SA;
			} else if (alg.compare("cve") == 0) {
				opt->algorithm = MERLIN_ALGO_MMAP_CVE;
			} else if (alg.compare("bn2cn") == 0) {
				opt->algorithm = MERLIN_ALGO_CONVERT;
			} else if (alg.compare("generator") == 0) {
				opt->algorithm = MERLIN_ALGO_GENERATOR;
			} else {
				std::string err_msg("Algorithm ");
				err_msg += alg + " is not supported.";
				throw std::runtime_error(err_msg);
			}
		}

		// parse the time limit
		if (vm.count("time-limit")) {
			opt->timeLimit = vm["time-limit"].as<int>();
		}

		// parse the random generator seed
		if (vm.count("seed")) {
			opt->seed = vm["seed"].as<size_t>();
		}

		// parse the number of iterations
		if (vm.count("iterations")) {
			opt->iterations = vm["iterations"].as<int>();
		}

		// parse the output format
		if (vm.count("output-format")) {
			std::string format = vm["output-format"].as<std::string>();
			if (format.compare("uai") == 0) {
				opt->outputFormat = MERLIN_OUTPUT_UAI;
			} else if (format.compare("json") == 0) {
				opt->outputFormat = MERLIN_OUTPUT_JSON;
			} else {
				std::string err_msg("The output format ");
				err_msg += format + " is not supported";
				throw std::runtime_error(err_msg);
			}
		}

		// parse epsilon value
		if (vm.count("epsilon")) {
			opt->epsilon = vm["epsilon"].as<double>();
		}

		// parse random flip probability
		if (vm.count("flip-proba")) {
			opt->flip_probability = vm["flip-proba"].as<double>();
		}

		// parse initialization method
		if (vm.count("init-method")) {
			opt->init_method = vm["init-method"].as<std::string>();
		}

		// parse initial temperature
		if (vm.count("init-temp")) {
			opt->init_temp = vm["init-temp"].as<double>();
		}

		// parse cooling factor
		if (vm.count("alpha")) {
			opt->alpha = vm["alpha"].as<double>();
		}

		// parse max number of flips
		if (vm.count("max-flips")) {
			opt->max_flips = vm["max-flips"].as<size_t>();
		}

		// parse the taboo list size
		if (vm.count("taboo-size")) {
			opt->taboo_size = vm["taboo-size"].as<size_t>();
		}

		// parse the cache size
		if (vm.count("cache-size")) {
			opt->cache_size = vm["cache-size"].as<size_t>();
		}

		// parse the number of nodes
		if (vm.count("nodes")) {
			opt->num_nodes = vm["nodes"].as<size_t>();
		}

		// parse the number of parents
		if (vm.count("parents")) {
			opt->num_parents = vm["parents"].as<size_t>();
		}

		// parse the number of instances
		if (vm.count("instances")) {
			opt->num_instances = vm["instances"].as<size_t>();
		}

		// parse the ksize
		if (vm.count("ksize")) {
			opt->ksize = vm["ksize"].as<size_t>();
		}

		// parse the percentage
		if (vm.count("kpercent")) {
			opt->kpercent = vm["kpercent"].as<size_t>();
		}

		// parse the graph type
		if (vm.count("graph-type")) {
			opt->graph_type = vm["graph-type"].as<std::string>();
		}

		// parse the graph type
		if (vm.count("query-type")) {
			opt->query_type = vm["query-type"].as<std::string>();
		}

		// parse the number of MAP (query) variables
		if (vm.count("num-query")) {
			opt->num_query = vm["num-query"].as<size_t>();
		}

		// parse the number of sample queries
		if (vm.count("num-samples")) {
			opt->num_samples = vm["num-samples"].as<size_t>();
		}

		// parse the number of extras (ie, extra LCN statements)
		if (vm.count("num-extras")) {
			opt->num_extras = vm["num-extras"].as<size_t>();
		}

		// parse the number of evidence (ie, LCN evidence)
		if (vm.count("num-evid")) {
			opt->num_evid = vm["num-evid"].as<size_t>();
		}

	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		delete opt;
		return NULL;
	}

	return opt;
}

