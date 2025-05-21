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

    // Initialize the query variables (MAP variables)
    m_query.clear();
    for (size_t v = 0; v < nvar(); ++v) {
        if (m_evidence.find(v) != m_evidence.end()) {
            continue;
        }

        m_query.push_back(v);
    }
}

// Build the weighted mini-bucket heuristic
double map2u::build_heuristic() {
    
    // Number of variables
    size_t num_vars = nvar();
    std::mt19937 rng(1234);
    num_vars += 1; // include the dummy
    m_buckets.clear();
    m_intermediate.clear();
    m_augmented.clear();
    
    std::cout << "[HEUR] Building the WMB heuristic..." << std::endl;

    // Initialize the buckets
    std::cout << "[HEUR] Initialize the buckets." << std::endl;
    std::vector<bool> used(num_vars, false);
    m_buckets.resize(num_vars);
    m_intermediate.resize(num_vars);
    m_augmented.resize(num_vars);
    for (size_t i = 0; i < m_order.size(); ++i) {
        size_t v = m_order[i];
        m_buckets[i].set_variable(v);
        for (size_t j = 0; j < m_factors.size(); ++j) {
            interval& f = m_factors[j];
            int ch = f.get_child();
            if (used[ch] == true) {
                continue;
            } else {
                // check if the current interval factor contains the bucket var
                if (f.vars().contains(var(v))) {
                    used[ch] = true;
                    // buckets[i].add_potential(f.to_potential(false));

                    potential p = f.to_potential(false);
                    p.approximate(m_potential_approx, m_potential_size, m_epsilon);
                    m_buckets[i].add_potential(p);
                }
            }
        }
    }

    // Eliminate the variables (following the elimination ordering)
    std::cout << "[HEUR] Begin variable elimination ..." << std::endl;
    bucket& scalars = m_buckets.back();
    for (size_t i = 0; i < num_vars - 1; ++i) {
        size_t v = m_order[i];
        variable vx = var(v);


        // Partition the bucket into mini-buckets
        std::vector<potential> partition = m_buckets[i].create_partition(
            m_ibound, m_potential_approx, m_potential_size, m_epsilon);

        if (m_verbose > 0) {
            std::cout << "[HEUR] Eliminating variable: " << v << std::endl;
            std::cout << "  - created " << partition.size() << " mini-buckets" << std::endl;
        }

        // Moment-matching between the mini-buckets
        if (m_matching_strategy > 0 && partition.size() > 1) { // match between multiple mini-buckets
            moment_matching(vx, partition);
        }

        // Eliminate the bucket variable from each mini-bucket
        for (size_t j = 0; j < partition.size(); ++j) {

            // Combine the potentials in the mini-bucket and eliminate the variable
            potential& result = partition[j];

            // Eliminate the variable (in-place)
            result.elim_max(vx);

            // Approximate the potential (in-place)
            result.approximate(m_potential_approx, m_potential_size, m_epsilon);

            // Remove dominated elements from the potential (max/min)
            if (m_query_type == MERLIN_MAP_MAXIMAX) {
                result.maximize();
            } else if (m_query_type == MERLIN_MAP_MAXIMIN) {
                result.minimize();
            }

            // Place new potential in the appropriate bucket
            if (result.nvar() == 0) { // i.e., scalar == empty scope
                scalars.add_potential(result);
                m_intermediate[num_vars - 1].push_back(result);
            } else {
                // Find the closest bucket that contains a variable in the potential's scope
                for (size_t j = i + 1; j < num_vars - 1; ++j) {
                    int y = m_buckets[j].get_variable();
                    variable vy = var(y);
                    if (result.vars().contains(vy)) {
                        m_buckets[j].add_potential(result);
                        m_augmented[y].push_back(result);
                        break;
                    } else {
                        m_intermediate[y].push_back(result);
                    }
                }
            }
        } // done mini-buckets

    } // done elimination

    std::cout << "[HEUR] Finished variable elimination." << std::endl;

    // After elimination, combine all scalars to determine global bound
    potential r(1.0);
    std::vector<potential>& pots = scalars.potentials();
    for (size_t i = 0; i < pots.size(); ++i) {
        r.multiply(pots[i]);
    }
    
    // Prune dominated scalars (no need for approximation -- just scalars)
    if (m_query_type == MERLIN_MAP_MAXIMAX) {
        r.maximize();
    } else if (m_query_type == MERLIN_MAP_MAXIMIN) {
        r.minimize();
    }

    // Check for singleton
    if (r.p().size() > 1) {
        std::cout << "[WMB] WARNING: more than one final scalars detected: " << r.p().size() << std::endl; 
    }

    // Get the best score
    double global_bound = r.p()[0][0];
    std::cout << "[HEUR] Global bound: " << global_bound << " (" << std::log10(global_bound) << ")" << std::endl;
    std::cout << "[HEUR] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl; 
    std::cout << "[HEUR] Finished building the heuristic." << std::endl;

    if (m_verbose > 0) {
        std::cout << "[DEBUG] Bucket structure:" << std::endl;
        for (size_t i = 0; i < m_buckets.size(); ++ i) {
            std::cout << "Bucket [" << m_buckets[i].get_variable() << "]" << std::endl;
            std::vector<potential>& pots = m_buckets[i].potentials(); 
            for (size_t j = 0; j < pots.size(); ++j) {
                std::cout << pots[j] << std::endl;
            } 
        }
        
        std::cout << "[DEBUG] Intermediate structure:" << std::endl;
        for (size_t i = 0; i < m_intermediate.size(); ++ i) {
            std::cout << "Intermediate [" << i << "]" << std::endl;
            std::vector<potential>& pots = m_intermediate[i]; 
            for (size_t j = 0; j < m_intermediate[i].size(); ++j) {
                std::cout << m_intermediate[i][j] << std::endl;
            } 
        }
    }

    return global_bound;
}

// Get the heuristic value for a variable during search given the current assignment
double map2u::get_heuristic(size_t var, std::map<size_t, size_t>& assignment, bool upper) {

	// variable 'var' is assumed to be already assigned (in 'assignment')
	double h = 1.0;

	// go over augmented and intermediate lists and combine all values
	for (size_t i = 0; i < m_augmented[var].size(); ++i) {
        h *= m_augmented[var][i].get_value(assignment, upper);
	}
	for (size_t i = 0; i < m_intermediate[var].size(); ++i) {
		h *= m_intermediate[var][i].get_value(assignment, upper);
	}

	return h;    
}

// Depth-First Search
void map2u::dfs() {

    // Prologue
    std::cout << "[DFS] Running Depth-First Search for MAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[DFS] Query type: maximax" << std::endl;
    } else {
        std::cout << "[DFS] Query type: maximin" << std::endl;
    }
    std::cout << "[DFS] Num MAP vars: " << m_query.size() << std::endl;
    std::cout << "[DFS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Depth-first search
    std::vector<int> best_config;
    double best_score = -1.0;
    bool timeout = false;
    size_t num_vars = m_query.size();
    size_t num_sols = 0, num_nodes = 0;

    // Enumerate all possible assignments of the MAP variables
    std::vector<int> values(num_vars, 0);
    values[num_vars - 1] = -1;
    int i;
    std::cout << "[DFS] Start search ...:" << std::endl;
    while (true) {

        // Enumerate "parent" variables.
        for (i = num_vars - 1; i >= 0; --i) {
            if (values[i] < 1) break;
            values[i] = 0;
        }

        if (i < 0) break;	// done;
        ++values[i];
        num_nodes++;

        // NOW: all MAP variables have a specific value combination.
        std::map<size_t, size_t> config;
        for (size_t j = 0; j < m_query.size(); ++j) {
            config[m_query[j]] = values[j];
        }

        // Evaluate the current MAP assignment
        double score = 1.0;
        for (std::vector<interval>::iterator ci = m_factors.begin(); ci != m_factors.end(); ++ci) {
            interval& f = *ci;
            interval::value v = f.get_value(config);
            if (m_query_type == MERLIN_MAP_MAXIMIN) {
                score *= v.first;
            } else {
                score *= v.second;
            }
        }

        if (score > best_score) {
            best_score = score;
            best_config = values;
            num_sols++;

            std::cout << "   - found better solution [" << best_score << " (" << std::log10(best_score) << ")]: ";
            std::copy(best_config.begin(), best_config.end(), std::ostream_iterator<int>(std::cout, " "));
            std::cout << std::endl;
        }

        if (m_verbose > 0) {
            std::cout << "SOL: [" << score << " (" << std::log10(score) << ")]: ";
            std::copy(values.begin(), values.end(), std::ostream_iterator<int>(std::cout, " "));
            std::cout << std::endl;
        }

        // Check for timeout
        double elapsed = (timeSystem() - m_start_time);
        if (m_time_limit > 0 && elapsed > m_time_limit) {
            std::cout << "  - TIMELIMT" << std::endl;
            timeout = true;
        }

        if (timeout) {
            break; // timout
        }
    }

    // Assemble the solution
    m_best_cost = best_score;
    m_best_config.resize(num_vars);
    for (size_t i = 0; i < best_config.size(); ++i) {
        m_best_config[i] = best_config[i];
    }

    std::cout << "[DFS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[DFS] Best cost: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[DFS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[DFS] Solutions found: " << num_sols << std::endl;
    std::cout << "[DFS] Number of nodes: " << num_nodes << std::endl;
    std::cout << "[DFS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_cost = best_score;
}

// Moment-matching (max) between the mini-buckets of a sigle bucket
void map2u::moment_matching(variable vx, std::vector<potential>& partition) {
    // Moment matching strategies:
    // 0 - no moment matching
    // 1 - single function (PLUB/PGLB)
    // 2 - exhaustive
    // 3 - ...

    if (m_matching_strategy == 0) { // no moment matching
        return; 
    } else if (m_matching_strategy == 1) { // single function
       
        // Do moment matching between the mini-buckets
        size_t R = partition.size();
        std::vector<factor> ftmp(R);         // compute geometric mean
        variable_set vs = partition[0].vars();  // on all mutual variables
        for (size_t i = 1; i < R; i++) {
            vs &= partition[i].vars();
        }

        // The auxiliary lambdas are the PLUB(1) approx of the max-marginals
        factor fmatch(vs, 1.0);
        for (size_t i = 0; i < R; i++) {
            potential marg = maxmarginal(partition[i], vs); // max-marginal on common vars
            if (m_query_type == MERLIN_MAP_MAXIMAX) {
                marg.plub(1);
            } else {
                marg.pglb(1);
            }
            ftmp[i] = marg[0]; // save the max-marginal of the mini-bucket
            fmatch *= ftmp[i];
        }

        fmatch ^= (1.0/R);         // and match each bucket to it
        for (size_t i = 0; i < R; i++) {
            factor f = (fmatch/ftmp[i]);
            potential pot(f);
            partition[i].multiply(pot);
        }

    } else if (m_matching_strategy == 2) {
        throw 2; // not implemented yet
    } else {
        throw 1; // not implemented yet
    }

}

// Credal Weighted Mini-Buckets for MAP (approximate)
void map2u::wmb() {

    // Initialize the solver
    std::cout << "[CWMB] Running Credal Weighted Mini-Buckets for MAP" << std::endl;
    if (m_query_type == MERLIN_MAP_MAXIMIN) {
        std::cout << "[CWMB] Query type: maximin MAP" << std::endl;
    } else {
        std::cout << "[CWMB] Query type: maximax MAP" << std::endl;
    }

    // Number of variables
    size_t num_vars = nvar();
    std::mt19937 rng(1234);
     
    // Create the minfill elimination ordering
    std::vector<size_t> elim_order;
    elim_order = order2();
    std::cout << "[CWMB] Elimination order: ";
    std::copy(elim_order.begin(), elim_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[CWMB] Induced width: " << m_width << std::endl;
    std::cout << "[CWMB] MB ibound: " << m_ibound << std::endl;
    std::cout << "[CWMB] Number of variables: " << num_vars << std::endl;
    std::cout << "[CQMB] Moment matching: " << m_matching_strategy << std::endl;

    // Initialize the buckets
    std::cout << "[CWMB] Initialize the buckets" << std::endl;
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
                    // buckets[i].add_potential(f.to_potential(false));

                    potential p = f.to_potential(false);
                    p.approximate(m_potential_approx, m_potential_size, m_epsilon);
                    buckets[i].add_potential(p);

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

    // Eliminate the variables (following the elimination ordering)
    std::cout << "[CWMB] Begin variable elimination ..." << std::endl;
    std::vector<potential> scalars;
    bool timeout = false;
    for (size_t i = 0; i < num_vars; ++i) {
        size_t v = elim_order[i];
        variable vx = var(v);
        std::string vtype = "MAX";
        std::cout << "[CWMB] Eliminating " << vtype << " variable: " << v << std::endl;

        // Partition the bucket into mini-buckets
        std::vector<potential> partition = buckets[i].create_partition(
            m_ibound, m_potential_approx, m_potential_size, m_epsilon);
        std::cout << "  - created " << partition.size() << " mini-buckets" << std::endl;

        // Moment-matching between the mini-buckets
        if (m_matching_strategy > 0 && partition.size() > 1) { // match between multiple mini-buckets
            if (m_verbose > 0) {
                std::cout << "[DEBUG] Partition before moment-matching:" << std::endl;
                for (size_t j = 0; j < partition.size(); ++j) {
                    std::cout << partition[j] << std::endl;
                }
            }
    
            moment_matching(vx, partition);

            if (m_verbose > 0) {
                std::cout << "[DEBUG] Partition after moment-matching:" << std::endl;
                for (size_t j = 0; j < partition.size(); ++j) {
                    std::cout << partition[j] << std::endl;
                }
            }    
        }


        // Eliminate the bucket variable from each mini-bucket
        for (size_t j = 0; j < partition.size(); ++j) {

            std::cout << "  - processing mini-bucket: " << j << std::endl;

            // Combine the potentials in the mini-bucket and eliminate the variable
            potential& result = partition[j];
            if (m_verbose > 0) {
                std::cout << "[DEBUG] Mini-bucket:" << std::endl;
                std::cout << result << std::endl;
            }

            // Eliminate the variable (in-place)
            result.elim_max(vx);

            if (m_verbose > 0) {
                std::cout << "[DEBUG] Result before pruning:" << std::endl;
                std::cout << result << std::endl;
            }

            // Approximate the potential
            result.approximate(m_potential_approx, m_potential_size, m_epsilon);

            // Remove dominated elements from the potential (max/min)
            if (m_query_type == MERLIN_MAP_MAXIMAX) {
                result.maximize();
            } else if (m_query_type == MERLIN_MAP_MAXIMIN) {
                result.minimize();
            }

            std::cout << "  - generated potential size: " << result.size() << std::endl;
            
            if (m_verbose > 0) {
                std::cout << "[DEBUG] Result after pruning:" << std::endl;
                std::cout << result << std::endl;
            }

            // Place new potential in the appropriate bucket
            if (result.nvar() == 0) { // i.e., scalar == empty scope
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

    std::cout << "[CWMB] Finished variable elimination." << std::endl;
    if (timeout) {
        std::cout << "[CWMB] Timeout: yes" << std::endl;
        return;
    }

    // After elimination, combine all scalars
    potential r(1.0);
    for (size_t i = 0; i < scalars.size(); ++i) {
        r.multiply(scalars[i]);
    }
    
    // Prune dominated scalars (no need for approximation -- just scalars)
    if (m_query_type == MERLIN_MAP_MAXIMAX) {
        r.maximize();
    } else if (m_query_type == MERLIN_MAP_MAXIMIN) {
        r.minimize();
    }

    // Check for singleton
    if (r.p().size() > 1) {
        std::cout << "[WMB] WARNING: more than one final scalars detected: " << r.p().size() << std::endl; 
    }

    // Get the best score
    m_best_cost = r.p()[0][0];
    if (m_verbose > 0) {
        std::cout << "Final constant potential is:" << std::endl << r << std::endl;
    }

    // /*
    std::cout << "[CWMB] Generating the MAP configuration (bottom-up) ..." << std::endl;
    // Compute the MAP assignment; going backwards in the ordering
    std::map<size_t, size_t> config;
    for (int i = num_vars - 1; i >= 0; --i) {
        size_t v = elim_order[i];

        std::cout << "[CWMB] Processing MAX variable: " << v << std::endl;
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
        std::cout << "[CWMB] Argmax for variable " << v << " is " << val << std::endl;

        // Check for timeout
        if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
            std::cout << "  - TIMELIMT" << std::endl;
            timeout = true;
            break;
        }
    }
    std::cout << "[CWMB] Finished generating the MAP configuration." << std::endl;
    // */
    if (!timeout) {
        // Assemble the solution
        m_best_config.resize(m_query.size());
        for (size_t i = 0; i < m_query.size(); ++i) {
            m_best_config[i] = config[m_query[i]];
            // m_best_config[i] = -1;
        }

        std::cout << "[CWMB] Best solution: ";
        std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "[CWMB] Best cost: " << m_best_cost << " (" << std::log10(m_best_cost) << ")" << std::endl;
        std::cout << "[CWMB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
        std::cout << "[CWMB] Timeout: no" << std::endl;
    } else {
        std::cout << "[CWMB] Timeout: yes" << std::endl;
    }
}

std::string map2u::to_string(variable_set &vars, std::map<size_t, size_t> &config) {
    std::stringstream ss;
    variable_set::const_iterator ci = vars.begin();
    for (; ci != vars.end(); ++ci) {
        size_t varx = *ci;
        size_t val = config.at(varx);
        ss << " " << varx << "=" << val; 
    }

    return ss.str();
}

void map2u::set_cache_context(search_node* n, const std::set<size_t>& ctxt) const {

    std::stringstream signature;
    for (std::set<size_t>::const_iterator si = ctxt.begin(); si != ctxt.end(); ++si) {
        signature << "x" << *si << "=" << m_assignment.at(*si) << ";";
    }

    std::string str_context = signature.str();
	n->set_context(str_context);
}

search_node* map2u::next_leaf() {

	search_node* node = next_node();
	while (node != NULL) {

		// check for time limit violation
		if (m_time_limit > 0 && timeSystem() - m_start_time > m_time_limit) {
			throw SEARCH_TIMEOUT;
		}

		if (do_process(node)) { // initial processing
			return node;
		}
		if (do_caching(node)) { // caching?
			return node;
		}
		if (do_pruning(node)) { // pruning?
			return node;
		}
		if (do_expand(node)) { // node expansion
			return node;
		}
		node = next_node();
	}

	return NULL;
}

search_node* map2u::next_node() {
	if (!m_stack.empty()) {
		search_node* n = m_stack.top();
		m_stack.pop();
        return n;
	} 
    
    return NULL;
}

bool map2u::do_process(search_node* n) {

    // Safety checks
	assert(n != NULL);
	if (n->get_type() == MERLIN_NODE_AND) {
		size_t var = n->get_variable();
		size_t val = n->get_value();
		m_assignment[var] = val; // record assignment

	} else { // NODE_OR
		// do nothing
	}

	return false; // default
}

// Retrieve an OR node from the cache if previously cached (context-based)
bool map2u::do_caching(search_node* n) {

    // Safety checks
	assert(n != NULL);
	int var = n->get_variable();
	pseudotree_node* ptnode = m_pseudotree->get_node(var);

	if (n->get_type() == MERLIN_NODE_AND) { // AND node -> reset associated adaptive cache tables

        // no caching applied

	} else { // OR node, try actual caching

        // No caching at root
		if (!ptnode->get_parent()) {
			return false;
        }

		if (ptnode->get_context().size() <= ptnode->get_parent()->get_context().size()) {

			// add cache context information
            set_cache_context(n, ptnode->get_context());

			// try to get value from cache
			try {
				// will throw int(UNKNOWN) if not found
				std::pair<double, std::vector<int>> entry = m_search_space->read(var, n->get_context());
				n->set_cost( entry.first ); // set value
				n->set_assignment( entry.second ); // set assignment
				n->set_leaf(true); // mark as leaf
				++m_cache_hits;

                if (m_verbose > 0) {
                    std::cout << "[CACHE] Found cached OR node: " << n->to_string() << std::endl;
                }

				return true;
			} catch (...) { // cache lookup failed
				n->set_cachable(); // mark for caching later
			}
		}
	} // if on node type

	return false; // default, no caching applied

}

bool map2u::do_pruning(search_node* n) {

    // Safety checks
	assert(n != NULL);

	if (can_prune(n)) {
		n->set_leaf(true);
        n->set_pruned();
        if (m_verbose > 0) {
            std::cout << "[PRUNE] Found pruned node: " << n->to_string() << std::endl;
        }

        if (n->get_type() == MERLIN_NODE_OR) {
			if (isnan(n->get_cost())) { // value could be set by LDS
				n->set_cost(0.0);
            }
		} else if (n->get_type() == MERLIN_NODE_AND) {
			n->set_cost(0.0); // dead end
		}

		return true;
	}

	return false; // default false
}

bool map2u::do_expand(search_node* n) {

    // Safety checks
	assert(n != NULL);
	std::vector<search_node*> expanded;

	if (n->get_type() == MERLIN_NODE_AND) {  // AND node

		// Update the heuristic
		std::map<size_t, size_t> assignment = n->get_path_assignment();
		size_t var = n->get_variable();
        // m_heuristic->update(var, assignment);

        // Generate the OR children of an AND node (if any)
		if (generate_children(n, expanded)) {
			return true; // no children
        }

        // Push children onto the stack
        std::vector<search_node*>::reverse_iterator it = expanded.rbegin();
		for (; it != expanded.rend(); ++it) {
			m_stack.push(*it);
        }

	} else if (n->get_type() == MERLIN_NODE_OR) {  // OR node

        // Generate the AND children of an OR node (if any)
		if (generate_children(n, expanded)) {
			return true; // no children
        }

        // Push children onto the stack
        std::vector<search_node*>::reverse_iterator it = expanded.rbegin();
		for (; it != expanded.rend(); ++it) {
			m_stack.push(*it);
		} // for loop

	} // if over node type

	return false; // default false (children generated)
}

double map2u::heuristic(search_node* n) {

    // Safety checks
    assert(n && n->get_type() == MERLIN_NODE_OR);

    // Get the OR node variable
	int var = n->get_variable();
	std::vector<double> dv;
    dv.resize(m_domains[var] * 2);

    bool upper = (m_query_type == MERLIN_MAP_MAXIMAX) ? true : false;
    double h = -INFINITY; // the new OR nodes h value
    std::map<size_t, size_t> assignment = m_assignment;
	std::list<potential>& funs = m_pseudotree->get_potentials(var);
	for (size_t k = 0; k < m_domains[var]; ++k) {
		assignment[var] = k;

		// compute heuristic value
		dv[2 * k] = get_heuristic(var, assignment, upper);

		// precompute weight value
		double w = 1.0;
        std::list<potential>::iterator li = funs.begin();
		for (; li != funs.end(); ++li) {
            potential& p = (*li);
			w *= p.get_value(assignment, upper);
		}

		// store label and heuristic into cache table
		dv[2 * k + 1] = w; // label
		dv[2 * k] *= w; // heuristic (includes label)

        if (dv[2 * k] > h) {
            h = dv[2 * h]; // keep max. for OR node heuristic (MAP var)
        }
	}

	n->set_heur(h);
	n->set_cache(dv);

	return h;    
}

bool map2u::generate_children(search_node* n, std::vector<search_node*>& chi) {
    
    // Safety checks
    assert(n != NULL);

    // Expand an AND node
    if (n->get_type() == MERLIN_NODE_AND) {

        // Get the AND node variable
        size_t var = n->get_variable();
        pseudotree_node* ptnode = m_pseudotree->get_node(var);

        // Increase AND node expansions
        m_search_space->add_node(MERLIN_NODE_AND);

        if (m_verbose > 0) {
            std::cout << "[EXPAND] Expanding AND node: " << n->to_string() << std::endl;
        }

        // Create new OR children (going in reverse due to reversal on stack)
        std::vector<pseudotree_node*>::const_reverse_iterator it = ptnode->get_children().rbegin();
        for (; it != ptnode->get_children().rend(); ++it) {

            // Get the pseudotree child
            int child_var = (*it)->get_variable();
            
            // Create the OR child
            search_node* c = new search_node(child_var, -1, MERLIN_NODE_OR);
            c->set_parent(n);
         
            // Compute and set heuristic estimate, includes child weights
            heuristic(c);
            c->set_depth(n->get_depth() + 1);
            chi.push_back(c);

            if (m_verbose > 0) {
                std::cout << "  - OR child: " << c->to_string() << std::endl;
            }
        } // for loop over new OR children

        if (chi.empty()) {
            n->set_leaf(true); // terminal node
            n->set_cost(1.0);
            return true; // no children
        }

        n->add_children(chi);

        return false; // default
    } else { // Expand an OR node
        assert(n && n->get_type() == MERLIN_NODE_OR);
       
        // Get the OR node variable
        int var = n->get_variable();
    
        // Increase OR node expansions
        m_search_space->add_node(MERLIN_NODE_OR);

        if (m_verbose > 0) {
            std::cout << "[EXPAND] Expanding OR node: " << n->to_string() << std::endl;
        }

        // Retrieve precomputed weights and heuristic values
        std::vector<double>& heur = n->get_cache();
        for (int val = m_domains[var] - 1; val >= 0; --val) {
            // early pruning if heuristic is zero (since it's an upper bound)
            if (heur[2 * val] == 0) { // 2*i=heuristic, 2*i+1=label
                continue;
            }
    
            search_node* c = new search_node(var, val, MERLIN_NODE_AND); // uses cached label
            c->set_parent(n);

            // Set cached heur. value (includes the weight)
            c->set_weight(heur[2 * val + 1]);
            c->set_heur(heur[2 * val]);
            c->set_depth(n->get_depth() + 1);
            
            chi.push_back(c);

            if (m_verbose) {
                std::cout << "  - AND child: "<< c->to_string() << std::endl;
            }
        }
    
        if (chi.empty()) { // deadend
            n->set_leaf(true);
            n->set_cost(0.0);
            return true; // no children
        }
    
        // sort new nodes by decreasing heuristic value - largest UB first
        // (use reverse iterator due to stack reversal)
        sort(chi.begin(), chi.end(), search_node::heur_greater);
    
        n->add_children(chi);
    
        return false; // default    
    }
} 

bool map2u::can_prune(search_node* n) {

    // Check if pruning is enabled
    if (!m_pruning) {
        return false; // disable pruning for now
    }

	// heuristic is an upper bound, hence can use to prune if value=0
	if (n->get_heur() == 0.0) {
		++m_num_deadends;
		return true;
	}

	search_node* curAND;
	search_node* curOR;
	double curPSTVal;

	if (n->get_type() == MERLIN_NODE_AND) {
		curAND = n;
		curOR = n->get_parent();
		curPSTVal = curAND->get_heur(); // includes label
	} else { // NODE_OR
		curAND = NULL;
		curOR = n;
		curPSTVal = curOR->get_heur(); // n->getHeur()
	}

	std::list<search_node*> notOptOR; // marks nodes for tagging as possibly not optimal

	// up to root node, if we have to
	while (curOR->get_parent()) {

		if ( curPSTVal <= curOR->get_cost() ) {
			for (std::list<search_node*>::iterator it = notOptOR.begin(); it != notOptOR.end(); ++it) {
				(*it)->set_optimal(false); // mark as possibly not optimal
            }

			++m_num_deadends;
			return true;// pruning is possible!
		}

		notOptOR.push_back(curOR);

		// climb up, update values
		curAND = curOR->get_parent();

		// collect AND node label
		curPSTVal *= curAND->get_weight();
		// incorporate already solved sibling OR nodes
		curPSTVal *= curAND->get_subsolved();
		// incorporate new not-yet-solved sibling OR nodes through their heuristic
		std::vector<search_node*>& children = curAND->get_children();
		for (size_t i = 0; i < children.size(); ++i) {
			if (!children[i] || children[i] == curOR) {
                continue;
            } else {
                curPSTVal *= children[i]->get_heur();
            }
		}
		curOR = curAND->get_parent();
	}

	// default, no pruning possible
	return false;
}

// Init the search space
search_node* map2u::init_search_space(double global_bound, double global_constant) {

    assert(m_search_space->get_root() == NULL);

	// Add initial set of dummy nodes.

	// create root OR node (dummy variable)
	pseudotree_node* ptroot = m_pseudotree->get_root();
    size_t root_var = ptroot->get_variable();
	search_node* root = new search_node(root_var, -1, MERLIN_NODE_OR);
	root->set_heur(global_bound);
	m_search_space->set_root(root);
    m_search_space->add_node(MERLIN_NODE_OR);

	// create dummy AND node (domain size 1) with global constant as label
	search_node* next = new search_node(root_var, 0, MERLIN_NODE_AND);
    m_search_space->add_node(MERLIN_NODE_AND);
    next->set_parent(root);
    next->set_weight(global_constant);
	root->add_child(next);
	next->set_heur( global_bound/next->get_weight() );

	return next;
}


/// Brute force search with exact CVE based evaluation (exact)
void map2u::bnb() {

    // Prologue
    std::cout << "[BB] Running OR Branch and Bound for Credal MAP" << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[BB] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[BB] Query type: maximin" << std::endl;
    } 

    std::cout << "[BB] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    // Number of variables
    size_t num_vars = nvar();
    std::mt19937 rng(1234);
    size_t dummy = num_vars; // dummy variable
    bool timeout = false;
    double* _EmergencyMem = new double[10]; // a memory buffer
    bool is_chain = (m_ao_search ? false : true); // chain pseudo tree (OR search)

    // Create the minfill elimination ordering (for precompiled heuristics)
    m_order = order2();
    std::cout << "[BB] Elimination order: ";
    std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[BB] Induced width: " << m_width << std::endl;
    std::cout << "[BB] MB ibound: " << m_ibound << std::endl;
    std::cout << "[BB] Number of variables: " << num_vars << std::endl;
    std::cout << "[BB] AND/OR search: " << (m_ao_search ? "yes" : "no") << std::endl;
    std::cout << "[BB] Enable caching: " << (m_caching ? "yes" : "no") << std::endl;
    std::cout << "[BB] Enable pruning: " << (m_pruning ? "yes" : "no") << std::endl;
    std::cout << "[BB] Moment matching: " << m_matching_strategy << std::endl;
    std::cout << "[BB] Chain PT: " << (is_chain ? "yes" : "no") << std::endl;

    // Moralize the graph
    graph g = this->moralize();

    // Branch and bound search
    m_solved = false;
    m_best_cost = -1;
    m_cache_hits = 0;
    m_num_deadends = 0;

    // Create the pseudo tree
    m_pseudotree = std::make_unique<pseudotree>();
    m_pseudotree->init(num_vars);
    m_pseudotree->build(g, m_order, is_chain);
    m_pseudotree->reset_potentials(m_factors); // includes the dummy variable ?
    m_order.push_back(dummy);

    // Output the pseudo tree
    if (m_verbose > 0) {
        m_pseudotree->dump(std::cout);
    }

    // Update the variable domains (including the dummy root of the pseudo tree)
    m_domains = m_dims;
    m_domains.push_back(1); // dummy var has 1 value!

    std::cout << "[BB] Pseudo tree width: " << m_pseudotree->get_width() << std::endl;
    std::cout << "[BB] Pseudo tree height: " << m_pseudotree->get_height() << std::endl;

    // Build the heuristic
    double global_bound = build_heuristic(); // get the global bound
    double global_constant = 1.0; 

    // Init search space
    m_search_space = std::make_unique<search_space>();
    m_search_space->init(num_vars);

    // Init the bound propagator (set caching as well)
    m_propagator = std::make_unique<bound_propagator>();
    m_propagator->init(m_start_time, m_pseudotree.get(), m_search_space.get(), m_caching);
    m_propagator->set_verbosity(m_verbose);

    std::cout << "[BB] Begin search..." << std::endl;

    try {

        // Init the search space
        search_node* first = init_search_space(global_bound, global_constant);
        if (first) {
            m_stack.push(first);
        }

        // Search
		search_node* n = next_leaf();
		while (n != NULL) { // throws timeout
			m_propagator->propagate(n, true); // true = report solutions
			m_best_cost = m_propagator->get_best_cost();
            m_best_config = m_propagator->get_best_config();
            n = next_leaf();
		}

	// 	// Proved optimality
		m_solved = true;
	} catch (std::bad_alloc& ba) {
		delete[] _EmergencyMem;
		_EmergencyMem = NULL;
		std::cout << "Critical out of memory exception! Aborting!";
	} catch (int e) {
        timeout = true;
    }

    size_t num_sols = m_propagator->get_num_solutions();
    std::pair<size_t, size_t> nodes_expanded = m_search_space->get_num_nodes();

    std::cout << "[BB] Finished search." << std::endl;
    std::cout << "[BB] Problem solved: " << (m_solved ? "true" : "false") << std::endl;
    std::cout << "[BB] Best solution: ";
    std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[BB] Best cost: " << m_best_cost << " (" << std::log10(m_best_cost) << ")" << std::endl;
    std::cout << "[BB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[BB] Number of AND nodes: " << nodes_expanded.first << std::endl;
    std::cout << "[BB] Number of OR nodes: " << nodes_expanded.second << std::endl;
    std::cout << "[BB] Cache hits: " << m_cache_hits << std::endl;
    std::cout << "[BB] Deadends: " << m_num_deadends << std::endl;
    std::cout << "[BB] Solutions found: " << num_sols << std::endl;
    std::cout << "[BB] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Clean up
    if (_EmergencyMem != NULL) {
        delete[] _EmergencyMem;
        _EmergencyMem = NULL;
    }
    if (m_search_space->get_root() != NULL) {
        delete m_search_space->get_root();
    }
}

std::vector<int> map2u::init_config() {
    std::vector<int> config;
    for (size_t i = 0; i < m_query.size(); ++i) {
        variable x = var(m_query[i]);
        //size_t val = randi2(x.states());
        size_t val = (randu() < 0.5 ? 0 : 1);
        config.push_back(val);
    }

    return config;
}

// Calculate the score of a MAP configuration: config is an assignment to all vars.
// It does not include the dummy variable.
double map2u::score(const std::vector<int>& config) {

    // Safety checks
    assert(config.size() == m_query.size());

    // NOW: all MAP variables have a specific value combination.
    std::map<size_t, size_t> assignment;
    for (size_t j = 0; j < m_query.size(); ++j) {
        assignment[j] = (size_t) config[j];
    }

    // Evaluate the current MAP assignment
    double val = 1.0;
    std::vector<potential>::iterator ci = m_potentials.begin();
    for (; ci != m_potentials.end(); ++ci) {
        potential& p = *ci;
        if (m_query_type == MERLIN_MAP_MAXIMAX) {
            double v = p.get_value(assignment, true); // upper
            val *= v;
        } else {
            double v = p.get_value(assignment, false); // lower
            val *= v;
        }
    }

    return val;
}

std::string map2u::make_key(const std::vector<int>& config) {
    std::ostringstream oss;
    for (size_t i = 0; i < m_query.size(); ++i) {
        oss << "x" << m_query[i] << "=" << config[i];
    }
    return oss.str();
}

void map2u::find_neighbors(const std::vector<int>& config, 
    std::vector<std::vector<int> >& neighbors) {

    neighbors.clear();
    for (size_t i = 0; i < m_query.size(); ++i) {
        variable x = var(m_query[i]);
        size_t num_states = x.states();
        for (size_t val = 0; val < num_states; ++val) {
            if (val != config[i]) {
                std::vector<int> new_config(config);
                new_config[i] = val;
                neighbors.push_back(new_config);
            }
        }
    }

    assert(neighbors.size() > 0);
}

// Stochastic Local Search
void map2u::sls() {
    // Init the cache
    std::map<std::string, double> cache;

    // Prologue
    std::cout << "[SLS] Running Stochastic Local Search for MAP" << std::endl;
    std::cout << "[SLS] Total iterations: " << m_iterations << std::endl;
    std::cout << "[SLS] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[SLS] Random flip probability: " << m_flip_probability << std::endl;
    std::cout << "[SLS] Max cached configs: " << m_cache_size << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[SLS] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[SLS] Query type: maximin" << std::endl;
    }
    std::cout << "[SLS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "[SLS] Potentials created: ";
    for (size_t j = 0; j < m_factors.size(); ++j) {
        interval& f = m_factors[j];
        potential p = f.to_potential(false);
        m_potentials.push_back(p);
    }
    std::cout << m_potentials.size() << std::endl;

    // Keep track of the overall best configuration
    std::vector<int> best_config, current_config;
    double best_score = -1.0, current_score = -1.0;
    size_t total_flips = 0, total_hits = 0;

    // Perform stochastic hill climbing for a number of iterations
    size_t num_sols = 0;
    bool timeout = false;
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

        std::cout << "[SLS] Iteration #" << iter << " ... " << std::endl;
        std::cout << "[SLS]   New initial solution: ";
        std::copy(current_config.begin(), current_config.end(), 
            std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "[SLS]   New initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;
        std::vector<int> best_config_iter;
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
            std::vector<int> next_config;
            double next_score = -1.0;

            // Neighbors of the current config
            std::vector<std::vector<int> > neighbors;
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
                std::vector<int> best_neighbor;
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

                std::cout << "[SLS]   - found better solution [" << best_score_iter << " (" << std::log10(best_score_iter) << ")" << "] after " << total_flips << " flips: ";
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

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "[SLS] TIMEOUT" << std::endl;
                timeout = true;
                break;
            }
        }
        
        // Check for timeout
        if (timeout) {
            break;
        }

        std::cout << "[SLS]   - finished after " << m_max_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[SLS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[SLS] Best cost: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[SLS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[SLS] Solutions found: " << num_sols << std::endl;
    std::cout << "[SLS] Total flips: " << total_flips << std::endl;
    std::cout << "[SLS] Total hits: " << total_hits << std::endl;
    std::cout << "[SLS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_cost = best_score;
}

// Taboo Search
void map2u::ts() {
    // Init the cache
    std::map<std::string, double> cache;

    // Initialize the taboo search
    std::cout << "[TS] Running Taboo Search for MAP" << std::endl;
    std::cout << "[TS] Total iterations: " << m_iterations << std::endl;
    std::cout << "[TS] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[TS] Taboo list size: " << m_taboo_size << std::endl;
    std::cout << "[TS] Max cached configs: " << m_cache_size << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[TS] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[TS] Query type: maximin" << std::endl;
    }
    std::cout << "[TS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "[TS] Potentials created: ";
    for (size_t j = 0; j < m_factors.size(); ++j) {
        interval& f = m_factors[j];
        potential p = f.to_potential(false);
        m_potentials.push_back(p);
    }
    std::cout << m_potentials.size() << std::endl;

    // Keep track of the overall best configuration
    std::vector<int> best_config, current_config;
    double best_score = -1.0, current_score = -1.0;

    // Create the taboo list
    std::map<std::string, bool> taboo_list; // keeps track of visited configs

    // Perform taboo search for a number of iterations
    size_t num_sols = 0, total_flips = 0, total_hits = 0;
    bool timeout = false;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Best config in the current iteration
        std::vector<int> best_config_iter;
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
        std::cout << "[TS] Iteration #" << iter << " ... " << std::endl;
        std::cout << "[TS]   New initial solution: ";
        std::copy(current_config.begin(), current_config.end(), 
            std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "[TS]   New initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;

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
            std::vector<int> best_neighbor;            
            std::vector<std::vector<int> > neighbors;
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
            std::vector<int> next_config;
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

                std::cout << "[TS]   - found better solution [" << best_score_iter << " (" << std::log10(best_score_iter) << ")" << "] after " << total_flips << " flips: ";
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

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "[TS] TIMEOUT" << std::endl;
                timeout = true;
                break;
            }
        }

        // Check for timeout
        if (timeout) {
            break;
        }

        std::cout << "[TS]   - finished after " << total_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[TS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[TS] Best cost: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[TS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[TS] Solutions found: " << num_sols << std::endl;
    std::cout << "[TS] Total flips: " << total_flips << std::endl;
    std::cout << "[TS] Total hits: " << total_hits << std::endl;
    std::cout << "[TS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_cost = best_score;

}

// Simulated Annealing
void map2u::sa() {
    // Init the cache
    std::map<std::string, double> cache;

    // Generate the initial configuration
    std::cout << "[SA] Running Simulated Annealing for MAP" << std::endl;
    std::cout << "[SA] Total iterations: " << m_iterations << std::endl;
    std::cout << "[SA] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[SA] Max cached configs: " << m_cache_size << std::endl;
    std::cout << "[SA] Initial temperature: " << m_init_temperature << std::endl;
    std::cout << "[SA] Cooling factor (alpha): " << m_alpha << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[SA] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[SA] Query type: maximin" << std::endl;
    }
    std::cout << "[SA] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "[TS] Potentials created: ";
    for (size_t j = 0; j < m_factors.size(); ++j) {
        interval& f = m_factors[j];
        potential p = f.to_potential(false);
        m_potentials.push_back(p);
    }
    std::cout << m_potentials.size() << std::endl;

    std::vector<int> current_config = init_config();
    double current_score = score(current_config);
    cache[make_key(current_config)] = current_score;

    std::cout << "[SA] Initial solution: ";
    std::copy(current_config.begin(), current_config.end(), 
        std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[SA] Initial score: " << current_score << " (" << std::log10(current_score) << ")" << std::endl;

    // Keep track of the overall best configuration
    std::vector<int> best_config = current_config;
    double best_score = current_score;
    size_t total_flips = 0, total_hits = 0;

    // Perform simulated annealing for a number of iterations (restart annealing)
    size_t num_sols = 0;
    bool timeout = false;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Restart annealing from the current best config
        std::vector<int> current_config = best_config;
        double current_score = best_score;
        double T = m_init_temperature;
        std::cout << "[SA] Iteration #" << iter << " ... " << std::endl;
        std::cout << "[SA]   - initial temperature: " << T << std::endl;
        // Perform simulated annealing for a max number of flips
        for (size_t flip = 1; flip <= m_max_flips; ++flip) {
            total_flips++;

            // Attempt to move to a random neighbor
            std::vector<int> next_config;
            double next_score = -1.0;
            std::vector<std::vector<int> > neighbors;
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

                std::cout << "[SA]   - found better solution [" << best_score << " (" << std::log10(best_score) << ")" << "] after " << flip << " flips and T " << T << ": ";
                std::copy(best_config.begin(), best_config.end(), 
                    std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
            }

            // Adjust the temperature
            // T *= m_alpha;
            if (flip % 100 == 0) {
                T *= m_alpha;
            }

            // Prune cache table if full
            while (cache.size() > m_cache_size) {
                cache.erase(cache.begin());
            }

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "[SA] TIMOUT" << std::endl;
                timeout = true;
                break;
            }
        }

        // Check for timeout
        if (timeout) {
            break;
        }

        std::cout << "[SA]   - final temperature: " << T << std::endl;
        std::cout << "[SA]   - finished after " << total_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[SA] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[SA] Best cost: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[SA] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[SA] Solutions found: " << num_sols << std::endl;
    std::cout << "[SA] Total flips: " << total_flips << std::endl;
    std::cout << "[SA] Total hits: " << total_hits << std::endl;
    std::cout << "[SA] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_cost = best_score;

}

void map2u::update_penalties(std::vector<int>& config) {

    // NOW: all MAP variables have a specific value combination.
    std::map<size_t, size_t> assignment;
    for (size_t j = 0; j < m_query.size(); ++j) {
        assignment[j] = (size_t) config[j];
    }

    // Find out the maximum utility
    double max_util = -infty();
    bool upper = (m_query_type == MERLIN_MAP_MAXIMAX ? true : false);
    for (size_t i = 0; i < m_potentials.size(); ++i) {
        potential& pot = m_potentials[i];
        factor& p = m_penalties[i];
        double util = -pot.get_value(assignment, upper) / (1.0 + p.get_value(assignment));
        max_util = std::max(max_util, util);
    }

    // Update the penalties
    for (size_t i = 0; i < m_penalties.size(); ++i) {
        potential& pot = m_potentials[i];
        factor& p = m_penalties[i];
        double util = -pot.get_value(assignment, upper) / (1.0 + p.get_value(assignment));
        if (util == max_util) {
            double v = p.get_value(assignment);
            p.set_value(assignment, v + 1.0);
        }
    }
}

void map2u::scale_penalties() {
    for (size_t i = 0; i < m_penalties.size(); ++i) {
        factor& p = m_penalties[i];
        p.scale(0.00001);
    }
}

// Calculate the score of a MAP configuration: config is an assignment to all vars.
// It does not include the dummy variable.
std::pair<double, double> map2u::score_gls(const std::vector<int>& config) {

    // Safety checks
    assert(config.size() == m_query.size());

    // NOW: all MAP variables have a specific value combination.
    std::map<size_t, size_t> assignment;
    for (size_t j = 0; j < m_query.size(); ++j) {
        assignment[j] = (size_t) config[j];
    }

    // Evaluate the probability of the current MAP assignment
    bool upper = (m_query_type == MERLIN_MAP_MAXIMAX ? true : false);
    double p_val = 1.0, g_val = 0.0, w = 1000.0;
    for (size_t i = 0; i < m_potentials.size(); ++i) {
        potential& pot = m_potentials[i];
        factor& p = m_penalties[i];
        double v = pot.get_value(assignment, upper);
        double l = p.get_value(assignment);
        double log_v = (v == 0.0 ? 10000.0 : std::log10(v));

        p_val *= v;
        g_val += (log_v - w*l);
    }

    return std::make_pair(p_val, g_val);
}

// Guided Local Search
void map2u::gls() {

    // Init the cache
    std::map<std::string, std::pair<double, double> > cache;

    // Prologue
    std::cout << "[GLS] Running Guided Local Search for MAP" << std::endl;
    std::cout << "[GLS] Total iterations: " << m_iterations << std::endl;
    std::cout << "[GLS] Flips per iteration: " << m_max_flips << std::endl;
    std::cout << "[GLS] Random flip probability: " << m_flip_probability << std::endl;
    std::cout << "[GLS] Max cached configs: " << m_cache_size << std::endl;
    if (m_query_type == MERLIN_MMAP_MAXIMAX) {
        std::cout << "[GLS] Query type: maximax" << std::endl;
    } else if (m_query_type == MERLIN_MMAP_MAXIMIN) {
        std::cout << "[GLS] Query type: maximin" << std::endl;
    }
    std::cout << "[GLS] Query vars: ";
    std::copy(m_query.begin(), m_query.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "[GLS] Create potentials: ";
    for (size_t j = 0; j < m_factors.size(); ++j) {
        interval& f = m_factors[j];
        potential p = f.to_potential(false);
        m_potentials.push_back(p);
    }
    std::cout << m_potentials.size() << std::endl;

    std::cout << "[GLS] Initialize penalties: ";
    for (size_t j = 0; j < m_factors.size(); ++j) {
        interval& f = m_factors[j];
        factor p(f.vars(), 0.0); // zero penalties
        m_penalties.push_back(p);
    }
    std::cout << m_penalties.size() << std::endl;

    // Keep track of the overall best configuration
    std::vector<int> best_config, current_config;
    double best_score = -infty(), current_score = -infty();
    double best_cost = -infty(), current_cost = -infty();
    size_t total_flips = 0, total_hits = 0;

    // Perform GLS for a number of iterations
    size_t num_sols = 0;
    bool timeout = false;
    for (size_t iter = 1; iter <= m_iterations; ++iter) {

        // Start with a new random initial config
        current_config = init_config();
        std::string ckey = make_key(current_config);
        std::map<std::string, std::pair<double, double> >::iterator mi = cache.find(ckey);
        if (mi != cache.end()) {
            total_hits++;
            std::pair<double, double> cval = mi->second;
            current_score = cval.second; // objective val
            current_cost = cval.first; // probability val
        } else {
            std::pair<double, double> cval = score_gls(current_config);
            current_score = cval.second; // objective val
            current_cost = cval.first; // probability val
            cache[ckey] = cval;
        }

        std::cout << "[GLS] Iteration #" << iter << " ... " << std::endl;
        std::cout << "[GLS]   New initial solution: ";
        std::copy(current_config.begin(), current_config.end(), 
            std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "[GLS]   New initial score: " << current_score << " [" << current_cost << " (" << std::log10(current_cost) << ")]" << std::endl;
        std::vector<int> best_config_iter;
        double best_score_iter = -infty(), best_cost_iter = -infty();

        // Keep track of the overall best configuration
        if (current_score > best_score) {
            best_config = current_config;
            best_score = current_score;
            best_cost = current_cost;
        }

        // Keep track of the best config during current iteration
        if (current_score > best_score_iter) {
            best_config_iter = current_config;
            best_score_iter = current_score;
            best_cost_iter = current_cost;
        }

        // Start flipping variables (one round of SLS)
        for (size_t flip = 1; flip <= m_max_flips; ++flip) {
            total_flips++;

            // Next config (neighbor) to move to
            std::vector<int> next_config;
            double next_score = -infty(), next_cost = -infty();

            // Neighbors of the current config
            std::vector<std::vector<int> > neighbors;
            find_neighbors(current_config, neighbors);

            // Toss the coin (p)
            double p = randu();
            if (p <= m_flip_probability) {
                // Select a random neighbor
                size_t j = randi2((int)neighbors.size()); // random neighbor
                next_config = neighbors[j];
                std::string ckey = make_key(next_config);
                std::map<std::string, std::pair<double, double> >::iterator ci = cache.find(ckey);
                if (ci != cache.end()) {
                    total_hits++;
                    std::pair<double, double> cval = ci->second;
                    next_score = cval.second;
                    next_cost = cval.first;
                } else {
                    std::pair<double, double> cval = score_gls(next_config);
                    next_score = cval.second;
                    next_cost = cval.first;
                    cache[ckey] = cval;
                }
            } else {
                // Find the best scoring neighbor
                double best_neighbor_score = -infty(), best_neighbor_cost = -infty();
                std::vector<int> best_neighbor;
                for (size_t k = 0; k < neighbors.size(); ++k) {
                    std::string nkey = make_key(neighbors[k]);
                    double neighbor_score = -1.0, neighbor_cost = -1.0;
                    std::map<std::string, std::pair<double, double> >::iterator ci = cache.find(nkey);
                    if (ci != cache.end()) {
                        std::pair<double, double> cval = ci->second;
                        neighbor_score = cval.second;
                        neighbor_cost = cval.first;
                    } else {
                        std::pair<double, double> cval = score_gls(neighbors[k]);
                        neighbor_score = cval.second;
                        neighbor_cost = cval.first;
                        cache[nkey] = cval;
                    }

                    if (neighbor_score > best_neighbor_score) {
                        best_neighbor_score = neighbor_score;
                        best_neighbor_cost = neighbor_cost;
                        best_neighbor = neighbors[k];
                    }
                }

                // If best scoring neighbor is better then keep it
                if (best_neighbor_score > current_score) {
                    next_config = best_neighbor;
                    next_score = best_neighbor_score;
                    next_cost = best_neighbor_cost;
                } else { // Otherwise, select a random neighbor
                    size_t j = randi2((int)neighbors.size()); // random neighbor
                    next_config = neighbors[j];
                    std::string ckey = make_key(next_config);
                    std::map<std::string, std::pair<double, double> >::iterator ci = cache.find(ckey);
                    if (ci != cache.end()) {
                        total_hits++;
                        std::pair<double, double> cval = ci->second;
                        next_score = cval.second;
                        next_cost = cval.first;
                    } else {
                        std::pair<double, double> cval = score_gls(next_config);
                        next_score = cval.second;
                        next_cost = cval.first;
                        cache[ckey] = cval;
                    }
                }
            }

            // Check if next config score is better in the current iteration
            if (next_score > best_score_iter) {
                best_score_iter = next_score;
                best_cost_iter = next_cost;
                best_config_iter = next_config;

                std::cout << "[SLS]   - found better solution: " << best_score_iter << " [" << best_cost_iter << " (" << std::log10(best_cost_iter) << ")]: ";
                std::copy(best_config_iter.begin(), best_config_iter.end(), 
                    std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
            } 

            // Keep track of the overall best config
            if (next_score > best_score) {
                best_score = next_score;
                best_cost = next_cost;
                best_config = next_config;
                num_sols++;                
            }

            // Move to the next config
            current_config = next_config;
            current_score = next_score;
            current_cost = next_cost;

            // Prune cache table if full
            while (cache.size() > m_cache_size) {
                cache.erase(cache.begin());
            }

            // Check for timeout
            if (m_time_limit > 0 && (timeSystem() - m_start_time) > m_time_limit) {
                std::cout << "[GLS] TIMEOUT" << std::endl;
                timeout = true;
                break;
            }
        } // end for
        
        // Update penalties
        update_penalties(best_config);

        // Scale penalties
        scale_penalties();

        // Check for timeout
        if (timeout) {
            break;
        }

        std::cout << "[GLS]   - finished after " << m_max_flips << " flips, " << total_hits << " hits and " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    }

    std::cout << "[GLS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[GLS] Best score: " << best_score << std::endl;
    std::cout << "[GLS] Best cost: " << best_cost << " (" << std::log10(best_cost) << ")" << std::endl;
    std::cout << "[GLS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[GLS] Solutions found: " << num_sols << std::endl;
    std::cout << "[GLS] Total flips: " << total_flips << std::endl;
    std::cout << "[GLS] Total hits: " << total_hits << std::endl;
    std::cout << "[GLS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_cost = best_score;

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
        bnb();
    } else if (m_search_method.compare("wmb") == 0) { // Weighted Mini-Buckets
        wmb();
    } else if (m_search_method.compare("sls") == 0) { // Stochastic Local Search
        sls();
    } else if (m_search_method.compare("sa") == 0) { // Simulated Annealing
        sa();
    } else if (m_search_method.compare("ts") == 0) { // Taboo Search
        ts();
    } else if (m_search_method.compare("gls") == 0) { // Guided Local Search
        gls();
    }
    
}

// Write the solution to the output stream
void map2u::write_solution(std::ostream& out, int output_format) {
	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
        out << " \"task\" : \"MAP\", ";
        out << " \"value\" : " << std::fixed
            << std::setprecision(MERLIN_PRECISION)
            << (m_best_cost) << ", ";
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
