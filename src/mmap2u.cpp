/*
 * mmap2u.cpp
 *
 *  Created on: 29 Sep 2021
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


#include "mmap2u.h"
#include "loopy2u.h"

namespace merlin {

// Initialize the solver
void mmap2u::init() {
	// Prologue
	std::cout << "[MMAP] Begin initialization ..." << std::endl;
	std::cout << "[MMAP] Random generator seed: " << m_seed << std::endl;
	rand_seed(m_seed); // set the random number generator seed

    // Init scorer for shc, ts, sa and gls
    if ( (m_search_method.compare("shc") == 0)
        || (m_search_method.compare("ts") == 0)
        || (m_search_method.compare("sa") == 0)
        || (m_search_method.compare("gls") == 0) ) {

        std::cout << "[MMAP] Initialize scorer ..." << std::endl;
        init_scorer();
    }
}

// Initialize scorer: apply HMM transformation to the MAP and evidence variables 
// by adding m auxiliary variables (where m is the number of MAP variables) and 
// m auxiliary deterministic factors in an HMM (poly-tree) structure.
void mmap2u::init_scorer() {
    
    // Create the parents set
    std::vector<size_t> parents(m_query);
    for (std::map<size_t, size_t>::iterator mi = m_evidence.begin(); 
            mi != m_evidence.end(); ++mi) {
        parents.push_back(mi->first);
    }

    // Create the auxiliary structure
    std::cout << "[MMAP] Create auxiliary structure for evidence ..." << std::endl;
    std::cout << "[MMAP] Number of auxiliary variables: " << parents.size() << std::endl;
    std::cout << "[MMAP] Evidence variables: ";
    std::copy(parents.begin(), parents.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Get the original factors
    std::vector<interval> fs = get_factors();
    
    // Add auxiliary variables and deterministic factors
    size_t x = nvar();
    std::vector<variable> prev;
    for (size_t i = 0; i < parents.size(); ++i) {
        variable ch = variable(x + i, 2);
        variable_set scope(ch);
        scope |= var(parents[i]);
        if (!prev.empty()) {
            scope |= prev.back();
        }

        interval f(scope);
        f.set_child(x + i);
        prev.push_back(ch);
        fs.push_back(f);

        m_aux_fid.push_back(fs.size() - 1);
        m_aux_vid.push_back(x + i);
    }

    if (m_verbose > 0) {
        std::cout << "[MMAP] Auxiliary variable indices: ";
        std::copy(m_aux_vid.begin(), m_aux_vid.end(), std::ostream_iterator<vindex>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "[MMAP] Auxiliary factor indices: ";
        std::copy(m_aux_fid.begin(), m_aux_fid.end(), std::ostream_iterator<findex>(std::cout, " "));
        std::cout << std::endl;
    }

    std::ostringstream oss;
    oss << "StopIter=" << m_iterations << ","
        << "Threshold=" << m_threshold << ","
        << "Verbose=" << 0
        << "Seed=" << m_seed;
    m_scorer = new loopy2u(fs);
    m_scorer->set_properties(oss.str());
}

// Update the auxiliary factors with deterministic information
void mmap2u::update_aux_factors(const std::map<size_t, size_t>& config) {

    for (size_t i = 0; i < m_aux_fid.size(); ++i) {
        findex fid = m_aux_fid[i];
        interval& f = m_scorer->get_factor(fid);
        interval::value val_0(0.0, 0.0), val_1(1.0, 1.0);
        index_config idx1(f.vars(), true);
        config_index idx2(f.vars(), true);
        vindex child = f.get_child();
        for (size_t j = 0; j < f.numel(); ++j) {
            std::map<size_t, size_t> asgn = idx1.convert(j);
            f[j] = (asgn[child] == 0) ? val_1 : val_0; // default
        }

        std::map<size_t, size_t> temp;
        for (variable_set::const_iterator it = f.vars().begin(); 
                it != f.vars().end(); ++it) {
            size_t a = it->label();
            if (a != child) {
                std::map<size_t, size_t>::const_iterator mi = config.find(a);
                if (mi != config.end()) {
                    temp[a] = mi->second;
                } else {
                    temp[a] = 1;
                }
            }
        }

        temp[child] = 1;
        f[idx2.convert(temp)] = val_1;
        temp[child] = 0;
        f[idx2.convert(temp)] = val_0;
    }

    if (m_verbose > 0) {
        std::cout << "[MMAP] Updated auxiliary factors:" << std::endl;
        for (size_t i = 0; i < m_aux_fid.size(); ++i) {
            const interval& f = m_scorer->get_factor(m_aux_fid[i]);
            std::cout << f << std::endl;
        }
    }
}

// Calculate the score of a MAP configuration. This is done by computing
// lower and upper bounds on P(e) where e is the evidence, namely, the current
// MAP assignment and the original evidence given in the input evidence file.
double mmap2u::score(const std::vector<size_t>& config) {

    // Safety checks
    assert(config.size() == m_query.size());

    if (m_verbose > 0) {
        std::cout << "[MMAP] Estimating score for config:";
        for (size_t j = 0; j < config.size(); ++j) {
            std::cout << " " << config[j];
        }
        std::cout << std::endl;
    }
    
    // Collect the original evidence and the MAP config
    variable_set parents;
    for (size_t i = 0; i < m_query.size(); ++i) {
        parents |= var(m_query[i]);
    }

    std::map<size_t, size_t> parents_config(m_evidence);
    for (size_t i = 0; i < m_query.size(); ++i) {
        size_t par = m_query[i];
        size_t val = config[i];
        parents_config[par] = val;
        parents |= var(par);
    }

    // To estimate the lower/upper bound on P(e) we extend the credal network by
    // adding a dummy variable X whose parents are the MAP variables and the
    // original evidence variables (if any). Then, we add a deterministic 
    // factor associated with X and pa(X) that is [1,1] iff pa(X) is consistent
    // with the MAP assignment and the original evidence. Subsequently, we run
    // Loopy2U (or IPE2U) to compute the marginal on X from and we only consider
    // the entry corresponding to X=1, namely the interval [l,u]: l <= P(e) <= u.
    // FIXME: this works only for up to 20-25 MAP variables; need another way
    //        to extend the CN so that we don't increase the induced width.
    update_aux_factors(parents_config);
    m_scorer->reset();
    m_scorer->run();

    vindex child = m_aux_vid.back();
    interval bel = m_scorer->belief(child);
    double lower = bel[1].first;
    double upper = bel[1].second;
    if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        return lower;
    } else if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        return upper;
    } else {
        return -1.0;
    }
}

// Compute the initial MAP configuration
std::vector<size_t> mmap2u::init_config() {
    if (m_init_method.compare("rand") == 0) {
        std::vector<size_t> config;
        for (size_t i = 0; i < m_query.size(); ++i) {
            variable x = var(m_query[i]);
            //size_t val = randi2(x.states());
            size_t val = (randu() < 0.5 ? 0 : 1);
            config.push_back(val);
        }

        return config;
    } else if (m_init_method.compare("mpe") == 0) {
        throw std::runtime_error("Initialization method MPE not implemented yet.");
    } else if (m_init_method.compare("mle") == 0) {
        throw std::runtime_error("Initialization method MLE not implemented yet.");
    } else {
        throw std::runtime_error("Unknown initialization method");
    }
}

// Get a random neighbor of the input MAP config
std::vector<size_t> mmap2u::get_random_neighbor(const std::vector<size_t>& config) {

    std::vector<size_t> neighbor(config);
    int i = randi2((int)m_query.size()); // select a MAP variable at random
    variable x = var(m_query[i]);
    size_t current_val = config[i];
    size_t next_val = 0;
    do {
        next_val = randi2(x.states()); // select a different value at random
    } while (next_val == current_val);

    assert (next_val != current_val);
    neighbor[i] = next_val;

    return neighbor;
}

// Find the neighbors of a MAP configuration conditioned on a variable
void mmap2u::find_neighbors(variable x, const std::vector<size_t>& config, 
    std::vector<std::vector<size_t> >& neighbors) {

    neighbors.clear();

    // Get the index of the variable
    std::vector<size_t>::const_iterator it = std::find(m_query.begin(), m_query.end(), x.label());
    assert(it != m_query.end());
    size_t i = it - m_query.begin();
    size_t num_states = x.states();
    for (size_t val = 0; val < num_states; ++val) {
        if (val != config[i]) {
            std::vector<size_t> new_config(config);
            new_config[i] = val;
            neighbors.push_back(new_config);
        }
    }

    assert(neighbors.size() > 0);
}

void mmap2u::find_neighbors(const std::vector<size_t>& config, 
    std::vector<std::vector<size_t> >& neighbors) {

    neighbors.clear();
    for (size_t i = 0; i < m_query.size(); ++i) {
        variable x = var(m_query[i]);
        size_t num_states = x.states();
        for (size_t val = 0; val < num_states; ++val) {
            if (val != config[i]) {
                std::vector<size_t> new_config(config);
                new_config[i] = val;
                neighbors.push_back(new_config);
            }
        }
    }

    assert(neighbors.size() > 0);
}

// Stochastic Hill Climbing (local search)
void mmap2u::hill_climbing2() {

    // Init the cache
    std::map<std::string, double> cache;

    // Prologue
    std::cout << "[HC] Running Stochastic Hill Climbing for MMAP" << std::endl;
    std::cout << "[HC] Initialization method: " << m_init_method << std::endl;
    std::cout << "[HC] Total iterations: " << m_iterations << std::endl;
    std::cout << "[HC] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[HC] Random flip probability: " << m_flip_probability << std::endl;
    std::cout << "[HC] Max cached configs: " << m_cache_size << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[HC] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[HC] Query type: maximin" << std::endl;
    } else {
        std::cout << "[HC] Query type: interval" << std::endl;
    }
    std::cout << "[HC] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Keep track of the overall best configuration
    std::vector<size_t> best_config, current_config;
    double best_score = -1.0, current_score = -1.0;
    size_t total_flips = 0, total_hits = 0;

    // Perform stochastic hill climbing for a number of iterations
    size_t num_sols = 0;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Start with a new random initial config
        current_config = init_config();
        std::string ckey = make_key(current_config);
        std::map<std::string, double>::iterator mi = cache.find(ckey);
        if (mi != cache.end()) {
            total_hits++;
            current_score = mi->second;
        } else {
            current_score = score(current_config);
            cache[ckey] = current_score;
        }

        std::cout << " Iteration #" << iter << " ... " << std::endl;
        std::cout << "   New initial solution: ";
        std::copy(current_config.begin(), current_config.end(), 
            std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "   New initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;
        std::vector<size_t> best_config_iter;
        double best_score_iter = -1.0;

        // Keep track of the overall best configuration
        if (current_score > best_score) {
            best_config = current_config;
            best_score = current_score;
        }

        // Keep track of the best config during current iteration
        if (current_score > best_score_iter) {
            best_config_iter = current_config;
            best_score_iter = current_score;
        }

        // Start flipping variables
        for (size_t flip = 1; flip <= m_max_flips; ++flip) {
            total_flips++;

            // Next config (neighbor) to move to
            std::vector<size_t> next_config;
            double next_score = -1.0;

            // Neighbors of the current config
            std::vector<std::vector<size_t> > neighbors;
            find_neighbors(current_config, neighbors);

            // Toss the coin (p)
            double p = randu();
            if (p <= m_flip_probability) {
                // Select a random neighbor
                size_t j = randi2((int)neighbors.size()); // random neighbor
                next_config = neighbors[j];
                std::string ckey = make_key(next_config);
                std::map<std::string, double>::iterator ci = cache.find(ckey);
                if (ci != cache.end()) {
                    total_hits++;
                    next_score = ci->second;
                } else {
                    next_score = score(next_config);
                    cache[ckey] = next_score;
                }
            } else {
                // Find the best scoring neighbor
                double best_neighbor_score = -1.0;
                std::vector<size_t> best_neighbor;
                for (size_t k = 0; k < neighbors.size(); ++k) {
                    std::string nkey = make_key(neighbors[k]);
                    double neighbor_score = -1.0;
                    std::map<std::string, double>::iterator ci = cache.find(nkey);
                    if (ci != cache.end()) {
                        neighbor_score = ci->second;
                    } else {
                        neighbor_score = score(neighbors[k]);
                        cache[nkey] = neighbor_score;
                    }

                    if (neighbor_score > best_neighbor_score) {
                        best_neighbor_score = neighbor_score;
                        best_neighbor = neighbors[k];
                    }
                }

                // If best scoring neighbor is better then keep it
                if (best_neighbor_score > current_score) {
                    next_config = best_neighbor;
                    next_score = best_neighbor_score;
                } else { // Otherwise, select a random neighbor
                    size_t j = randi2((int)neighbors.size()); // random neighbor
                    next_config = neighbors[j];
                    std::string ckey = make_key(next_config);
                    std::map<std::string, double>::iterator ci = cache.find(ckey);
                    if (ci != cache.end()) {
                        total_hits++;
                        next_score = ci->second;
                    } else {
                        next_score = score(next_config);
                        cache[ckey] = next_score;
                    }
                }
            }

            // Check if next config score is better in the current iteration
            if (next_score > best_score_iter) {
                best_score_iter = next_score;
                best_config_iter = next_config;

                std::cout << "   - found better solution [" << best_score_iter << " (" << std::log10(best_score_iter) << ")" << "] after " << total_flips << " flips: ";
                std::copy(best_config_iter.begin(), best_config_iter.end(), 
                    std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
            } 

            // Keep track of the overall best config
            if (next_score > best_score) {
                best_score = next_score;
                best_config = next_config;
                num_sols++;                
            }

            // Move to the next config
            current_config = next_config;
            current_score = next_score;

            // Prune cache table if full
            while (cache.size() > m_cache_size) {
                cache.erase(cache.begin());
            }
        }
        std::cout << "   - finished after " << m_max_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[HC] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[HC] Best score: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[HC] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[HC] Solutions found: " << num_sols << std::endl;
    std::cout << "[HC] Total flips: " << total_flips << std::endl;
    std::cout << "[HC] Total hits: " << total_hits << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_score = best_score;
}

std::string mmap2u::make_key(const std::vector<size_t>& config) {
    std::ostringstream oss;
    for (size_t i = 0; i < m_query.size(); ++i) {
        oss << "x" << m_query[i] << "=" << config[i];
    }
    return oss.str();
}

// Convert a configuration to a string
std::string mmap2u::to_string(variable_set& vars, std::map<size_t, size_t>& config) {
    std::ostringstream oss;
    for (variable_set::const_iterator ci = vars.begin(); ci != vars.end(); ++ci) {
        size_t var = ci->label();
        oss << "x" << var << "=" << config[var];
    }

    return oss.str();
}

// Tabu Search (local search)
void mmap2u::taboo_search2() {

    // Init the cache
    std::map<std::string, double> cache;

    // Initialize the taboo search
    std::cout << "[TS] Running Taboo Search for MMAP" << std::endl;
    std::cout << "[TS] Initialization method: " << m_init_method << std::endl;
    std::cout << "[TS] Total iterations: " << m_iterations << std::endl;
    std::cout << "[TS] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[TS] Taboo list size: " << m_taboo_size << std::endl;
    std::cout << "[TS] Max cached configs: " << m_cache_size << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[TS] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[TS] Query type: maximin" << std::endl;
    } else {
        std::cout << "[TS] Query type: interval" << std::endl;
    }
    std::cout << "[TS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Keep track of the overall best configuration
    std::vector<size_t> best_config, current_config;
    double best_score = -1.0, current_score = -1.0;

    // Create the taboo list
    std::map<std::string, bool> taboo_list; // keeps track of visited configs

    // Perform taboo search for a number of iterations
    size_t num_sols = 0, total_flips = 0, total_hits = 0;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Best config in the current iteration
        std::vector<size_t> best_config_iter;
        double best_score_iter = -1.0;

        // Generate a new random intial configuration
        current_config = init_config();
        taboo_list.clear();
        std::string ckey = make_key(current_config);
        std::map<std::string, double>::iterator mi = cache.find(ckey);
        if (mi != cache.end()) {
            current_score = mi->second;
            total_hits++;
        } else {
            current_score = score(current_config);
            cache[ckey] = current_score;
        }
        std::cout << " Iteration #" << iter << " ... " << std::endl;
        std::cout << "   New initial solution: ";
        std::copy(current_config.begin(), current_config.end(), 
            std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "   New initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;

        // Keep track of the overall best configuration
        if (current_score > best_score) {
            best_config = current_config;
            best_score = current_score;
        }

        // Keep track of the best config in the current iteration
        if (current_score > best_score_iter) {
            best_config_iter = current_config;
            best_score_iter = current_score;
        }

        // Repeat for a number of max flips per iteration
        for (size_t flip = 1; flip <= m_max_flips; ++flip) {

            total_flips++;

            // Add current config to taboo list (if enough space)
            assert(taboo_list.size() <= m_taboo_size);            
            std::string key = make_key(current_config);
            taboo_list[key] = true;

            // Find the best neighbor NOT in the taboo list
            double best_neighbor_score = -1.0;
            std::vector<size_t> best_neighbor;            
            std::vector<std::vector<size_t> > neighbors;
            find_neighbors(current_config, neighbors);
            for (size_t k = 0; k < neighbors.size(); ++k) {
                double neighbor_score = -1.0;
                std::string nkey = make_key(neighbors[k]);
                
                // Check if neighbor in taboo list (ignore if yes)
                std::map<std::string, bool>::iterator ti = taboo_list.find(nkey);
                if (ti == taboo_list.end()) {
                    // Compute the score
                    std::map<std::string, double>::iterator ci = cache.find(nkey);
                    if (ci != cache.end()) {
                        total_hits++;
                        neighbor_score = ci->second;
                    } else {
                        neighbor_score = score(neighbors[k]);
                        cache[nkey] = neighbor_score;
                    }

                    // Keep track of the best scoring neighbor not in taboo list
                    if (neighbor_score > best_neighbor_score) {
                        best_neighbor_score = neighbor_score;
                        best_neighbor = neighbors[k];
                    }
                }
            }

            // If no such neighbor exists then select one at random
            std::vector<size_t> next_config;
            double next_score = -1.0;
            if (best_neighbor_score == -1.0) {
                size_t j = randi2((int)neighbors.size()); // random neighbor
                next_config = neighbors[j];
                std::string ckey = make_key(next_config);
                std::map<std::string, double>::iterator ci = cache.find(ckey);
                if (ci != cache.end()) {
                    next_score = ci->second;
                    total_hits++;
                } else {
                    next_score = score(next_config);
                    cache[ckey] = next_score;
                }
            } else { // Found the best neighbor not in taboo list
                next_config = best_neighbor;
                next_score = best_neighbor_score;
            }

            // Keep track of best config in current iteration
            if (next_score > best_score_iter) {
                best_score_iter = next_score;
                best_config_iter = next_config;

                std::cout << "   - found better solution [" << best_score_iter << " (" << std::log10(best_score_iter) << ")" << "] after " << total_flips << " flips: ";
                std::copy(best_config.begin(), best_config.end(), 
                    std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
            }

            // Keep track of the overall best config
            if (next_score > best_score) {
                best_config = next_config;
                best_score = next_score;
                num_sols++;
            }

            // Move to the next config
            current_config = next_config;
            current_score = next_score;

            // Prune taboo list if full
            if (taboo_list.size() > m_taboo_size) {
                taboo_list.erase(taboo_list.begin());
            }

            // Prune cache table if full
            while (cache.size() > m_cache_size) {
                cache.erase(cache.begin());
            }
        }

        std::cout << "   - finished after " << total_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[TS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[TS] Best score: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[TS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[TS] Solutions found: " << num_sols << std::endl;
    std::cout << "[TS] Total flips: " << total_flips << std::endl;
    std::cout << "[TS] Total hits: " << total_hits << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_score = best_score;
}

// Simulated Annealing (local search)
void mmap2u::simulated_annealing2() {

    // Init the cache
    std::map<std::string, double> cache;

    // Generate the initial configuration
    std::cout << "[SA] Running Simulated Annealing for MMAP" << std::endl;
    std::cout << "[SA] Initialization method: " << m_init_method << std::endl;
    std::cout << "[SA] Total iterations: " << m_iterations << std::endl;
    std::cout << "[SA] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[SA] Max cached configs: " << m_cache_size << std::endl;
    std::cout << "[SA] Initial temperature: " << m_init_temperature << std::endl;
    std::cout << "[SA] Cooling factor (alpha): " << m_alpha << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[SA] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[SA] Query type: maximin" << std::endl;
    } else {
        std::cout << "[SA] Query type: interval" << std::endl;
    }
    std::cout << "[SA] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    std::vector<size_t> current_config = init_config();
    double current_score = score(current_config);
    cache[make_key(current_config)] = current_score;

    std::cout << "[SA] Initial solution: ";
    std::copy(current_config.begin(), current_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[SA] Initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;

    // Keep track of the overall best configuration
    std::vector<size_t> best_config = current_config;
    double best_score = current_score;
    size_t total_flips = 0, total_hits = 0;

    // Perform simulated annealing for a number of iterations (restart annealing)
    size_t num_sols = 0;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Restart annealing from the current best config
        std::vector<size_t> current_config = best_config;
        double current_score = best_score;
        double T = m_init_temperature;
        std::cout << " Iteration #" << iter << " ... " << std::endl;
        std::cout << "   - initial temperature: " << T << std::endl;
        // Perform simulated annealing for a max number of flips
        for (size_t flip = 1; flip <= m_max_flips; ++flip) {
            total_flips++;

            // Attempt to move to a random neighbor
            std::vector<size_t> next_config;
            double next_score = -1.0;
            std::vector<std::vector<size_t> > neighbors;
            find_neighbors(current_config, neighbors);
            size_t j = randi2((int)neighbors.size()); // random neighbor
            next_config = neighbors[j];
            std::string ckey = make_key(next_config);
            std::map<std::string, double>::iterator ci = cache.find(ckey);
            if (ci != cache.end()) {
                next_score = ci->second;
                total_hits++;
            } else {
                next_score = score(next_config);
                cache[ckey] = next_score;
            }

            // Compute Metropolis acceptance criterion (mac)
            double delta = std::log10(next_score) - std::log10(current_score);            
            if (delta > 0) { // next config is better; accept it
                current_config = next_config;
                current_score = next_score;
            } else {
                double p = randu();
                double threshold = std::exp(delta/T); // Metropolis acceptance criterion
                if (p < threshold) {
                    current_config = next_config; // move to a worse config
                    current_score = next_score;
                }
            }

            // Keep track of best config in the current iteration
            if (current_score > best_score) {
                best_config = current_config;
                best_score = current_score;
                num_sols++;

                std::cout << "   - found better solution [" << best_score << " (" << std::log10(best_score) << ")" << "] after " << flip << " flips: ";
                std::copy(best_config.begin(), best_config.end(), 
                    std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
            }

            // Adjust the temperature
            if (flip % 100 == 0) {
                T *= m_alpha;
            }

            // Prune cache table if full
            while (cache.size() > m_cache_size) {
                cache.erase(cache.begin());
            }
        }

        std::cout << "   - final temperature: " << T << std::endl;
        std::cout << "   - finished after " << total_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[SA] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[SA] Best score: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[SA] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[SA] Solutions found: " << num_sols << std::endl;
    std::cout << "[SA] Total flips: " << total_flips << std::endl;
    std::cout << "[SA] Total hits: " << total_hits << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_score = best_score;
}

// Guided Local Search
void mmap2u::guided_local_search2() {

}

// Variable Elimination for MMAP (exact)
void mmap2u::variable_elimination2() {

    // Initialize the solver
    std::cout << "[CVE] Running Credal Variable Elimination for MMAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[CVE] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[CVE] Query type: maximin" << std::endl;
    } else {
        std::cout << "[CVE] Interval query is not supported" << std::endl;
        std::cout << "[CVE] Stop" << std::endl;
        return;
    }
    std::cout << "[CVE] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Number of variables
    size_t num_vars = nvar();

    // Create constrained minfill ordering
    std::vector<size_t> elim_order;
    elim_order = constrained_order2(m_query);
    std::cout << "[CVE] Elimination order: ";
    std::copy(elim_order.begin(), elim_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[CVE] Induced width: " << m_width << std::endl;

    // Initialize the variable types
    std::vector<bool> var_types(num_vars, false);
    for (size_t i = 0; i < m_query.size(); ++i) {
        var_types[m_query[i]] = true; // mark as MAP variable
    }

    // Initialize the buckets
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
    for (size_t i = 0; i < num_vars; ++i) {
        size_t v = elim_order[i];
        variable vx = var(v);
        std::string vtype = (var_types[v] ? "MAX" : "SUM");
        std::cout << "[CVE] Eliminating " << vtype << " variable: " << v << std::endl;

        // Combine the potentials in the bucket
        potential comb(1.0);
        std::vector<potential>& pots = buckets[i].potentials();
        for (size_t j = 0; j < pots.size(); ++j) {
            comb.multiply(pots[j]);
        }

        // Eliminate the bucket variable
        potential result = comb;
        if (var_types[v] == true) { // MAX variable
            result.max(vx);
        } else {
            result.sum(vx);
        }

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
    }

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
        std::cout << "[CVE] WARNING: more than one final scalars detected: " << r.p().size() << std::endl; 
    }

    // Get the best score
    m_best_score = r.p()[0][0];

    // Compute the MAP assignment
    std::map<size_t, size_t> config;
    for (size_t i = num_vars - 1; i >= 0; --i) {
        size_t v = elim_order[i];
        if (var_types[v] == false) {
            break; // stop at the first SUM variable
        }

        std::cout << "[CVE] Processing MAX variable: " << v << std::endl;
        variable vx = var(v);
        potential result(1.0);
        std::vector<potential>& pots = buckets[i].potentials();
        for (size_t j = 0; j < pots.size(); ++j) {
            potential temp = pots[j];
            temp.substitute(config);
            result.multiply(temp);
        }

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Combined potential (before pruning):" << std::endl;
            std::cout << result << std::endl;
        }

        if (m_query_type == MERLIN_MMAP_MAXIMAX) {
            result.maximize();
        } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
            result.minimize();
        }

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Combined potential (after pruning):" << std::endl;
            std::cout << result << std::endl;

        }

        size_t val = result.argmax();
        config[v] = val;
        std::cout << "[CVE] Argmax for variable " << v << " is " << val << std::endl;
    }

    // Assemble the solution
    m_best_config.resize(m_query.size());
    for (size_t i = 0; i < m_query.size(); ++i) {
        m_best_config[i] = config[m_query[i]];
    }

    std::cout << "[CVE] Best solution: ";
    std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[CVE] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
    std::cout << "[CVE] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
}

// Run solver
void mmap2u::run() {

    // Init the start time
    m_start_time = timeSystem();

	// Initialize the algorithm
	init();

    // Run the search algorithms
    if (m_search_method.compare("hc") == 0) { // hill climbing
        hill_climbing2();
    } else if (m_search_method.compare("ts") == 0) { // tabu search
        taboo_search2();
    } else if (m_search_method.compare("sa") == 0) { // simulated annealing
        simulated_annealing2();
    } else if (m_search_method.compare("cve") == 0) { // credal variable elimination 
        variable_elimination2();
    }
}

// Write the solution to the output stream
void mmap2u::write_solution(std::ostream& out, int output_format) {
	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
        out << " \"task\" : \"MMAP\", ";
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
        out << "MMAP" << std::endl;
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
