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


#include "map2u.h"
#include "loopy2u.h"
#include "cve2u.h"

namespace merlin {

// Initialize the solver
void map2u::init() {
	// Prologue
	std::cout << "[MAP] Begin initialization ..." << std::endl;
	std::cout << "[MAP] Random generator seed: " << m_seed << std::endl;
	rand_seed(m_seed); // set the random number generator seed

    // Init scorer for shc, ts, sa and gls
    if ( (m_search_method.compare("bnb") == 0)
        || (m_search_method.compare("aobb") == 0)
        || (m_search_method.compare("bfs") == 0)
        || (m_search_method.compare("aobf") == 0) ) {

        std::cout << "[MAP] Precompile weighted mini-bucket heuristics ..." << std::endl;
        precompile_heuristics();
    }
}

void map2u::precompile_heuristics() {

}

// Depth-First Search
void map2u::dfs() {

    // Init the cache
    std::map<std::string, double> cache;

    // Prologue
    std::cout << "[DFS] Running Depth-First Search for MAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[DFS] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[DFS] Query type: maximin" << std::endl;
    } else {
        std::cout << "[DFS] Query type: interval" << std::endl;
    }
    std::cout << "[DFS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Keep track of the overall best configuration
    std::vector<size_t> best_config, current_config;
    double best_score = -1.0, current_score = -1.0;
    size_t total_flips = 0, total_hits = 0;
    size_t num_sols = 0;
    bool timeout = false;


    std::cout << "[DFS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[DFS] Best score: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[DFS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[DFS] Solutions found: " << num_sols << std::endl;
    std::cout << "[DFS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_score = best_score;
}

// Credal Weighted Mini-Buckets for MAP (approximate)
void map2u::wmb() {

    // Initialize the solver
    std::cout << "[WMB] Running Credal Weighted Mini-Buckets for MAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[WMB] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[WMB] Query type: maximin" << std::endl;
    } else {
        std::cout << "[WMB] Interval query is not supported" << std::endl;
        std::cout << "[WMB] Stop" << std::endl;
        return;
    }
    std::cout << "[CWMB] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Number of variables
    size_t num_vars = nvar();

    // Create constrained minfill ordering
    std::vector<size_t> elim_order;
    elim_order = constrained_order2(m_query);
    std::cout << "[WMB] Elimination order: ";
    std::copy(elim_order.begin(), elim_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[WMB] Induced width: " << m_width << std::endl;
    std::cout << "[WMB] MB ibound: " << m_ibound << std::endl;

    // Initialize the buckets
    std::cout << "[WMB] Initialize the buckets" << std::endl;
    std::vector<bool> used(num_vars, false);
    std::vector<bucket> buckets(num_vars);
    for (size_t i = 0; i < elim_order.size(); ++i) {
        size_t v = elim_order[i];
        buckets[i].set_variable(v);
        for (size_t j = 0; j < m_factors.size(); ++j) {
            interval& f = m_factors[j];
            int ch = f.get_child();
            if (used[ch] == true) {
                continue;
            } else {
                // check if the current interval factor contains the bucket var
                if (f.vars().contains(var(v))) {
                    used[ch] = true;
                    buckets[i].add_potential(f.to_potential());
                }
            }
        }
    }

    if (m_verbose > 0) {
        std::cout << "[DEBUG] Bucket structure:" << std::endl;
        for (size_t i = 0; i < buckets.size(); ++ i) {
            std::cout << "Bucket [" << buckets[i].get_variable() << "]" << std::endl;
            std::vector<potential>& pots = buckets[i].potentials(); 
            for (size_t j = 0; j < pots.size(); ++j) {
                std::cout << pots[j] << std::endl;
            } 
        }
    }

    // Eliminate the variables
    std::vector<potential> scalars;
    bool timeout = false;
    for (size_t i = 0; i < num_vars; ++i) {
        size_t v = elim_order[i];
        variable vx = var(v);
        std::string vtype = "MAX";
        std::cout << "[WMB] Eliminating " << vtype << " variable: " << v << std::endl;

        // Partition the bucket into mini-buckets
        std::vector<potential> partition = buckets[i].create_partition(m_ibound);
        std::cout << "  - created " << partition.size() << " mini-buckets" << std::endl;
        bool first = true;
        for (size_t j = 0; j < partition.size(); ++j) {

            std::cout << "  - processing mini-bucket: " << j << std::endl;

            // Combine the potentials in the mini-bucket and eliminate the variable
            potential& result = partition[j];
            if (m_verbose > 0) {
                std::cout << "[DEBUG] Mini-bucket:" << std::endl;
                std::cout << result << std::endl;
            }

            result.max(vx);

            if (m_verbose > 0) {
                std::cout << "[DEBUG] Result before pruning:" << std::endl;
                std::cout << result << std::endl;
            }

            // Remove dominated vertices
            if (m_query_type == MERLIN_MMAP_MAXIMAX) {
                result.maximize();
            } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
                result.minimize();
            }

            std::cout << "  - generated potential size: " << result.size() << std::endl;
            
            if (m_verbose > 0) {
                std::cout << "[DEBUG] Result after pruning:" << std::endl;
                std::cout << result << std::endl;
            }

            // Place new potential in the appropriate bucket
            if (result.isscalar()) {
                scalars.push_back(result);
            } else {
                // Find the closest bucket that contains a variable in the potential's scope
                for (size_t j = i + 1; j < num_vars; ++j) {
                    int y = buckets[j].get_variable();
                    variable vy = var(y);
                    if (result.vars().contains(vy)) {
                        buckets[j].add_potential(result);
                        break;
                    }
                }
            }

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "  - TIMELIMT" << std::endl;
                timeout = true;
                break;
            }
        } // done mini-buckets

        if (timeout) {
            break;
        }
    } // done elimination

    if (!timeout) {
        // After elimination, combine all scalars
        potential r(1.0);
        for (size_t i = 0; i < scalars.size(); ++i) {
            r.multiply(scalars[i]);
        }
        
        // Prune dominated scalars
        if (m_query_type == MERLIN_MMAP_MAXIMAX) {
            r.maximize();
        } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
            r.minimize();
        }

        // Check for singleton
        if (r.p().size() > 1) {
            std::cout << "[WMB] WARNING: more than one final scalars detected: " << r.p().size() << std::endl; 
        }

        // Get the best score
        m_best_score = r.p()[0][0];

        // Compute the MAP assignment
        std::map<size_t, size_t> config;
        for (size_t i = num_vars - 1; i >= 0; --i) {
            size_t v = elim_order[i];

            std::cout << "[WMB] Processing MAX variable: " << v << std::endl;
            variable vx = var(v);
            potential result(1.0);
            std::vector<potential>& pots = buckets[i].potentials();
            std::cout << "  - potentials in bucket: " << pots.size() << std::endl; 
            for (size_t j = 0; j < pots.size(); ++j) {
                potential temp = pots[j];
                if (m_verbose > 0) {
                    std::cout << "Before substitution:" << std::endl;
                    std::cout << temp << std::endl;
                }
                temp.substitute(config);
                if (m_verbose > 0) {
                    std::cout << "After substitiution:" << std::endl;
                    std::cout << temp << std::endl;
                }
                result.multiply(temp);
            }

            if (m_verbose > 0) {
                std::cout << "[DEBUG] Combined potential (before pruning):" << std::endl;
                std::cout << result << std::endl;
            }

            size_t val = result.argmax();
            config[v] = val;
            std::cout << "[WMB] Argmax for variable " << v << " is " << val << std::endl;

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "  - TIMELIMT" << std::endl;
                timeout = true;
                break;
            }
        }

        if (!timeout) {
            // Assemble the solution
            m_best_config.resize(m_query.size());
            for (size_t i = 0; i < m_query.size(); ++i) {
                m_best_config[i] = config[m_query[i]];
            }

            std::cout << "[WMB] Best solution: ";
            std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
            std::cout << std::endl;
            std::cout << "[WMB] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
            std::cout << "[WMB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
            std::cout << "[WMB] Timeout: no" << std::endl;
        } else {
            std::cout << "[WMB] Timeout: yes" << std::endl;
        }
    } else {
        std::cout << "[WMB] Timeout: yes" << std::endl;
    }
}

/// Brute force search with exact CVE based evaluation (exact)
void map2u::bnb() {

    // Prologue
    std::cout << "[BNB] Running Branch and Bound for MAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[BNB] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[BNB] Query type: maximin" << std::endl;
    } else {
        std::cout << "[BNB] Query type: interval" << std::endl;
    }
    std::cout << "[BNB] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Initialize the solver
    size_t num_vars = m_query.size();
    size_t num_sols = 0;
    double best_score = -1;
    std::vector<int> best_config;

    // Initialize the exact scorer (cve2u)
    std::ostringstream oss;
    oss << "Verbose=0,Seed=" << m_seed;
    std::vector<interval> fs = get_factors();
    merlin::cve2u exact_scorer(fs);
    exact_scorer.set_properties(oss.str());
    exact_scorer.init();

    // Enumerate all possible assignments of the MAP variables
    std::vector<int> values(num_vars, 0);
    values[num_vars - 1] = -1;
    int i;
    bool timeout = false;
    std::cout << "[BNB] Start search ...:" << std::endl;
    while (true) {

        // Enumerate "parent" variables.
        for (i = num_vars - 1; i >= 0; --i) {
            if (values[i] < 1) break;
            values[i] = 0;
        }

        if (i < 0) break;	// done;
        ++values[i];

        // NOW: all guery variables have a specific value combination.
        std::map<size_t, size_t> config;
        for (size_t j = 0; j < m_query.size(); ++j) {
            config[m_query[j]] = values[j];
        }

        // Evaluate the current MAP assignment
        std::pair<double, double> result = exact_scorer.eval(config);
        if (m_query_type == MERLIN_MMAP_MAXIMAX) {
            if (result.second > best_score) {
                best_score = result.second;
                best_config = values;
                num_sols++;

                std::cout << "   - found better solution [" << best_score << " (" << std::log10(best_score) << ")]: ";
                std::copy(best_config.begin(), best_config.end(), std::ostream_iterator<int>(std::cout, " "));
                std::cout << std::endl;
            }
        } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
            if (result.first > best_score) {
                best_score = result.first;
                best_config = values;
                num_sols++;

                std::cout << "   - found better solution [" << best_score << " (" << std::log10(best_score) << ")]: ";
                std::copy(best_config.begin(), best_config.end(), std::ostream_iterator<int>(std::cout, " "));
                std::cout << std::endl;
            }
        } else {
            // do nothing for now, but later collect the non-dominated ones
        }

        // Check for timeout
        double elapsed = (timeSystem() - m_start_time);
        if (m_time_limit > 0 && elapsed > m_time_limit) {
            std::cout << "  - TIMELIMT" << std::endl;
            timeout = true;
        }
    }

    // Assemble the solution
    m_best_score = best_score;
    m_best_config.resize(num_vars);
    for (size_t i = 0; i < best_config.size(); ++i) {
        m_best_config[i] = best_config[i];
    }

    std::cout << "[BNB] Finished search" << std::endl;
    std::cout << "[BNB] Best solution: ";
    std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[BNB] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
    std::cout << "[BNB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[BNB] Solutions found: " << num_sols << std::endl;
    std::cout << "[BNB] Timeout: " << (timeout ? "yes" : "no") << std::endl;
}

// Run solver
void map2u::run() {

    // Init the start time
    m_start_time = timeSystem();

	// Initialize the algorithm
	init();

    // Run the search algorithms
    if (m_search_method.compare("dfs") == 0) { // Depth-First Search
        dfs();
    } else if (m_search_method.compare("bnb") == 0) { // Branch and Bound Search
        bnb();
    } else if (m_search_method.compare("aobb") == 0) { // AND/OR Branch and Bound Search
        aobb();
    }
}

// Write the solution to the output stream
void map2u::write_solution(std::ostream& out, int output_format) {
	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
        out << " \"task\" : \"MAP\", ";
        out << " \"value\" : " << std::fixed
            << std::setprecision(MERLIN_PRECISION)
            << (m_best_score) << ", ";
        out << " \"status\" : \"true\", ";
        out << " \"solution\" : [ ";

        // Evidence variables are a disjoint set from the query variables
        for (vindex i = 0; i < m_query.size(); ++i) {
            vindex j = m_query[i];
            out << "{";
            out << " \"variable\" : " << j << ",";
            out << " \"value\" : " << m_best_config[i];
            out << "}";
            if (i != m_query.size() - 1) {
                out << ", ";
            }
        }
        out << "]}\n";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
        // evidence variables are a disjoint set from the query variables
        out << "MAP" << std::endl;
        out << m_query.size();
        for (vindex i = 0; i < m_query.size(); ++i) {
            vindex j = m_query[i];
            out << " " << j << " " << m_best_config[i];
        }
        out << std::endl;
	} else {
		std::string err_msg("[ERROR] Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}


} // end namespace
