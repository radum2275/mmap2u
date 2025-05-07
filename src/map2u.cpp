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

// AND/OR Branch and Bound search
void map2u::aobb() {

    bool timeout = false;
    size_t num_sols = 0;

    std::cout << "[AOBB] Running AND/OR Branch and Bound search ..." << std::endl;

    std::cout << "[AOBB] Finished search" << std::endl;
    std::cout << "[AOBB] Best solution: ";
    std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[AOBB] Best cost: " << m_best_cost << " (" << std::log10(m_best_cost) << ")" << std::endl;
    std::cout << "[AOBB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[AOBB] Solutions found: " << num_sols << std::endl;
    std::cout << "[AOBB] Timeout: " << (timeout ? "yes" : "no") << std::endl;
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
