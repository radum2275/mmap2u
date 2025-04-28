/*
 * merlin.cpp
 *
 *  Created on: 20 August 2018
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

#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "merlin.h"
#include "program_options.h"
#include "quickhull.h"

std::string fileToString(std::string filename) {
	std::ostringstream oss(std::ios::out | std::ios::binary); // *** binary
	std::ifstream in(filename.c_str());

	std::string line;
	while (std::getline(in, line)) {
		oss << line << std::endl;
	}

	return oss.str();
}

int test_quickhull() {
	using F = float;
	using value_type = double;
    constexpr std::size_t dim = 5;
    using Points = std::vector<std::array<F, dim>>;

    Points points(1500); // input

    { // fill it somehow (use real data)
        std::mt19937 gen;
        for (auto & [x, y, z, t, u] : points) {
            x = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
            y = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
            z = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
			t = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
			u = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
        }
    }

    const auto eps = std::numeric_limits<F>::epsilon();
    quick_hull<typename Points::const_iterator> qh{dim, eps};
    qh.add_points(std::cbegin(points), std::cend(points));
    auto initial_simplex = qh.get_affine_basis();
    if (initial_simplex.size() < dim + 1) {
		std::cout << "Degenerated input set." << std::endl;
        return 1; // degenerated input set
    }
    qh.create_initial_simplex(std::cbegin(initial_simplex), std::prev(std::cend(initial_simplex)));
    qh.create_convex_hull();
    if (!qh.check()) {
		std::cout << "Resulted structure is not convex (generally due to precision errors)." << std::endl;
        return 1; // resulted structure is not convex (generally due to precision errors)
    }

    // qh.facets_; // use as result
	std::cout << "Finished QuickHull." << std::endl;
	std::cout << "QH size (facets): " << qh.facets_.size() << std::endl;
    for (auto const & facet_ : qh.facets_) {
        auto const & vertices_ = facet_.vertices_;
        for (auto const & vertex_ : vertices_) {
            for (value_type const & coordinate_ : *vertex_) {
                //std::cout << coordinate_ << ' ';
            }
            //std::cout << '\n';
        }
	}

	return 0;
}

int main(int argc, char** argv) {

	ProgramOptions* opt = parseCommandLine(argc, argv);
	if (!opt) {
		std::cerr << "Invalid command line arguments." << std::endl;
		return MERLIN_EXIT_FAILURE;
	}

	// Set the output format if not set
	if (opt->outputFormat < 0) {
		opt->outputFormat = MERLIN_OUTPUT_UAI;
	}

	// Setup Merlin engine
	Merlin eng;
	eng.set_output_format(opt->outputFormat);
	eng.set_model_file(opt->modelFile);
	eng.set_evidence_file(opt->evidenceFile);
	eng.set_output_file(opt->outputFile);
	eng.set_query_file(opt->queryFile);
	eng.set_task(opt->task);
	eng.set_algorithm(opt->algorithm);
	eng.set_time_limit(opt->timeLimit);
	eng.set_scorer(opt->scorer);
	eng.set_ibound(opt->ibound);
	eng.set_iterations(opt->iterations);
	eng.set_debug(opt->debug);
	eng.set_verbose(opt->verbose);
	eng.set_threshold(opt->threshold);
	eng.set_epsilon(opt->epsilon);
	eng.set_seed(opt->seed);
	eng.set_flip_probability(opt->flip_probability);
	eng.set_init_method(opt->init_method);
	eng.set_init_temp(opt->init_temp);
	eng.set_alpha(opt->alpha);
	eng.set_max_flips(opt->max_flips);
	eng.set_taboo_size(opt->taboo_size);
	eng.set_cache_size(opt->cache_size);
	eng.set_num_nodes(opt->num_nodes);
	eng.set_num_parents(opt->num_parents);
	eng.set_num_instances(opt->num_instances);
	eng.set_ksize(opt->ksize);
	eng.set_kpercent(opt->kpercent);
	eng.set_graph_type(opt->graph_type);
	eng.set_query_type(opt->query_type);
	eng.set_num_query(opt->num_query);
	eng.set_num_samples(opt->num_samples);
	eng.set_num_extras(opt->num_extras);
	eng.set_num_evid(opt->num_evid);
	eng.set_potential_approx(opt->potential_approx);
	eng.set_potential_size(opt->potential_size);
	eng.set_moment_matching(opt->moment_matching);
	eng.set_caching(opt->caching);
	eng.set_pruning(opt->pruning);
	eng.set_force_bivalued(opt->force_bivalued);
	
	// Run the inference
	eng.init();
	int status = eng.run();

	delete opt;

	// test_quickhull();

	return status;
}

