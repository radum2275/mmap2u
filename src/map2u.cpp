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

    // Precompile heuristics
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

        // NOW: all MAP variables have a specific value combination.
        num_nodes++;
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

    std::cout << "[DFS] Best solution: ";
    std::copy(best_config.begin(), best_config.end(), 
        std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl << "[DFS] Best score: " << best_score << " (" << std::log10(best_score) << ")" << std::endl;
    std::cout << "[DFS] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[DFS] Solutions found: " << num_sols << std::endl;
    std::cout << "[DFS] Number of nodes: " << num_nodes << std::endl;
    std::cout << "[DFS] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Save best solution (and score)
    m_best_config = best_config;
    m_best_score = best_score;
}

// Moment-matching (max) in a mini-buckets partition
void map2u::moment_matching(variable vx, std::vector<potential>& partition, std::mt19937& rng) {

    // Special case: all mini-buckets are size 1 potentials (e.g., LPUB/GPLB(1))
    bool all_ones = true;
    for (size_t i = 0; i < partition.size(); ++i) {
        if (partition[i].size() > 1) {
            all_ones = false;
            break;
        }
    }

    size_t num_iters = m_iterations;
    if (all_ones) {
        num_iters = 1;
    }

    // At each iterations, select randomly functions from each potential
    for (size_t iter = 1; iter <= num_iters; ++iter) {
        std::vector<size_t> selection;
        for (size_t i = 0; i < partition.size(); ++i) {
            size_t sel = partition[i].select(rng);
            selection.push_back(sel);
        }

        // Do moment matching between the functions in the current selection
        size_t R = selection.size();
        std::vector<factor> ftmp(R);       // compute geometric mean
        variable_set vs = partition[0][selection[0]].vars(); // on all mutual variables
        for (size_t i = 1; i < R; i++) {
            size_t j = selection[i];
            vs &= partition[i][j].vars();
        }

        factor fmatch(vs, 1.0);
        for (size_t i = 0; i < R; i++) {
            size_t j = selection[i];
            ftmp[i] = maxmarginal(partition[i][j], vs);
            fmatch *= ftmp[i];
        }

        fmatch ^= (1.0/R);         // and match each bucket to it
        for (size_t i = 0; i < R; i++) {
            size_t j = selection[i];
            partition[i][j] *= (fmatch/ftmp[i]);
        }
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
                    p.approximate(m_query_type, m_potential_approx, m_potential_size, m_epsilon);
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
        std::vector<potential> partition = buckets[i].create_partition(m_ibound, m_query_type, m_potential_approx, m_potential_size, m_epsilon);
        std::cout << "  - created " << partition.size() << " mini-buckets" << std::endl;

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Partition before moment-matching:" << std::endl;
            for (size_t j = 0; j < partition.size(); ++j) {
                std::cout << partition[j] << std::endl;
            }
        }

        // Moment-matching between the mini-buckets
        if (m_matching && partition.size() > 1) { // match between multiple mini-buckets
            moment_matching(vx, partition, rng);
        }

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Partition after moment-matching:" << std::endl;
            for (size_t j = 0; j < partition.size(); ++j) {
                std::cout << partition[j] << std::endl;
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

            // Remove dominated vertices or approximate the potential
            result.approximate(m_query_type, m_potential_approx, m_potential_size, m_epsilon);

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
    m_best_score = r.p()[0][0];

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
        std::cout << "[CWMB] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
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

std::unique_ptr<search_node> map2u::next_leaf() {

	std::unique_ptr<search_node> node = std::move(next_node());
	while (node != nullptr) {

		// check for time limit violation
		if (m_time_limit > 0 && timeSystem() - m_start_time > m_time_limit) {
			throw SEARCH_TIMEOUT;
		}

		if (do_process(node.get())) { // initial processing
			return node;
		}
		if (do_caching(node.get())) { // caching?
			return node;
		}
		if (do_pruning(node.get())) { // pruning?
			return node;
		}
		if (do_expand(node.get())) { // node expansion
			return node;
		}
		node = std::move(next_node());
	}

	return nullptr;
}

std::unique_ptr<search_node> map2u::next_node() {
	if (m_stack.size() > 0) {
		std::unique_ptr<search_node> n = std::move(m_stack.top());
		m_stack.pop();
        return n;
	} 
    
    return nullptr;
}

bool map2u::do_process(search_node* n) {
	assert(n != NULL);
	if (n->get_type() == NODE_AND) {
		size_t var = n->get_variable();
		size_t val = n->get_value();
		m_assignment[var] = val; // record assignment

	} else { // NODE_OR
		// do nothing
	}

	return false; // default
}

bool map2u::do_caching(search_node* n) {
    return false;
}

bool map2u::do_pruning(search_node* n) {
    return false; 
}

bool map2u::do_expand(search_node* n) {
	assert(n);
	std::vector<std::unique_ptr<search_node>> expanded;

	if (n->get_type() == NODE_AND) {  // AND node

		// Update the heuristic
		std::vector<int> assignment;
		assignment.resize(m_assignment.size(), -1);
		n->get_path_assignment(assignment);
		size_t curr_var = n->get_variable();
        m_heuristic->update(assignment);

        // Generate the OR children of an AND node (if any)
		if (generate_children(n, expanded)) {
			return true; // no children
        }

        std::vector<std::unique_ptr<search_node>>::reverse_iterator it = expanded.rbegin();
		for (; it != expanded.rend(); ++it) {
			m_stack.push(std::move(*it));
        }

	} else if (n->get_type() == NODE_0R) {  // OR node

        // Generate the AND children of an OR node (if any)
		if (generate_children(n, expanded)) {
			return true; // no children
        }

        std::vector<std::unique_ptr<search_node>>::reverse_iterator it = expanded.rbegin();
		for (; it != expanded.rend(); ++it) {
			m_stack.push(std::move(*it));
		} // for loop

	} // if over node type

	return false; // default false    
}

// DONE: radu
bool map2u::generate_children(search_node* n, std::vector<std::unique_ptr<search_node>>& chi) {
    assert(n != nullptr);

    // Expand an AND node
    if (n->get_type() == NODE_AND) {
        assert(n && n->get_type() == NODE_AND);

        size_t var = n->get_variable();
        pseudotree_node* ptnode = m_pseudotree->get_node(var);

        // Increase AND node expansions
        m_num_nodes.first += 1;

        // Create new OR children (going in reverse due to reversal on stack)
        std::vector<pseudotree_node*>::const_reverse_iterator it = ptnode->get_children().rbegin();
        for (; it != ptnode->get_children().rend(); ++it) {

            // Get the pseudotree child
            int vChild = (*it)->get_variable();
            
            // Create the OR child
            std::unique_ptr<search_node> c = std::make_unique<search_node>(vChild, -1, NODE_OR);
          
            // Compute and set heuristic estimate, includes child labels
            heuristic(c);
            c->set_depth(n->get_depth() + 1);
            chi.push_back(std::move(c));

        } // for loop over new OR children

        if (chi.empty()) {
            n->set_leaf(); // terminal node
            n->set_cost(1);
            return true; // no children
        }

        // order subproblems in decreasing order of their heuristic - largest UB first
        // (use reverse iterator due to stack reversal)
        std::sort(chi.begin(), chi.end(), search_node::heur_greater);

        n->add_cildren(chi);

        return false; // default
    } else { // Expand an OR node
        assert(n->get_type() == NODE_OR);
       
        int var = n->get_variable();
    
        // Increase OR node expansions
        m_num_nodes.second += 1;
    
        // retrieve precomputed labels and heuristic values
        std::vector<double>& heur = n->get_heur_cache();
        for (int val = m_domains[var] - 1; val >= 0; --val) {
            // early pruning if heuristic is zero (since it's an upper bound)
            if (heur[2 * var] == 0) { // 2*i=heuristic, 2*i+1=label
                continue;
            }
    
            std::unique_ptr<search_node> c = std::make_unique<search_node>(var, val, NODE_AND); // uses cached label
            // set cached heur. value (includes the weight)
            c->set_weight(heur[2 * val + 1]);
            c->set_heur(heur[2 * val]);
            c->set_depth(n->get_depth() + 1);
            chi.push_back(c);
        }
    
        if (chi.empty()) { // deadend
            n->set_leaf();
            n->set_cost(0);
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
    return false;
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
    size_t num_sols = 0;
    double best_score = -1;
    bool timeout = false;
    std::vector<int> best_config(num_vars, -1);
    double* _EmergencyMem = new double[10]; // a memory buffer

    // Create the minfill elimination ordering (for precompiled heuristics)
    std::vector<size_t> elim_order;
    elim_order = order2();
    std::cout << "[BB] Elimination order: ";
    std::copy(elim_order.begin(), elim_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[BB] Induced width: " << m_width << std::endl;
    std::cout << "[BB] MB ibound: " << m_ibound << std::endl;
    std::cout << "[BB] Number of variables: " << num_vars << std::endl;

    // Branch and bound search
    m_solved = false;
	m_num_nodes = std::make_pair(0, 0);
    m_pseudotree = std::make_unique<pseudotree>();

    try {

		// Init the root of the search space
 		m_stack.push(std::make_unique<search_node>(num_vars, 0, NODE_AND));
        m_propagator = std::make_unique<bound_propagator>();

        // Search
		std::unique_ptr<search_node> n = next_leaf();
		while (n != nullptr) { // throws timeout
			m_propagator->propagate(n.get(), true); // true = report solutions
			m_best_score = m_propagator->get_score();
			n = next_leaf();
		}

		// Proved optimality
		m_solved = true;
	} catch (std::bad_alloc& ba) {
		delete[] _EmergencyMem;
		_EmergencyMem = NULL;
		std::cout << "Critical out of memory exception! Aborting!";
	} catch (int e) {
        timeout = true;
    }

    // Assemble the solution
    m_best_score = best_score;
    m_best_config.resize(num_vars);
    for (size_t i = 0; i < best_config.size(); ++i) {
        m_best_config[i] = best_config[i];
    }

    std::cout << "[BB] Finished search" << std::endl;
    std::cout << "[BB] Best solution: ";
    std::copy(m_best_config.begin(), m_best_config.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "[BB] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
    std::cout << "[BB] CPU time: " << (timeSystem() - m_start_time) << " seconds" << std::endl;
    std::cout << "[BB] Solutions found: " << num_sols << std::endl;
    std::cout << "[BB] Timeout: " << (timeout ? "yes" : "no") << std::endl;

    // Clean up
    if (_EmergencyMem != NULL) {
        delete[] _EmergencyMem;
        _EmergencyMem = NULL;
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
    std::cout << "[AOBB] Best score: " << m_best_score << " (" << std::log10(m_best_score) << ")" << std::endl;
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
        aobb();
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
