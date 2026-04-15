/*
 * bound_propagator.h
 *
 *  Created on: 17 Mar 2025
 *      Author: radu
 *
 * Copyright (c) 2025, International Business Machines Corporation. All rights reserved.
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


#include "bound_propagator.h"
#include "utils.h" 
namespace merlin {
 
// Constructor
bound_propagator::bound_propagator() {
    m_caching = false; // default no caching
    m_verbose = 0; // verbosity level
}

bound_propagator::~bound_propagator() {

}

/**
 * AND nodes do multiplication, while OR nodes do maximization.
 */

 search_node* bound_propagator::propagate(search_node* n, bool solution, search_node* upper_limit) {

    // 'n' is a leaf node (which triggers the propagation)
    assert(n->num_children() == 0);

    // These two pointers move upward in the search space, always one level
    // apart s.t. cur is the parent node of prev
    search_node* cur = n->get_parent(), *prev = n;
	if (cur == NULL) { // n is actually the root
        m_best_cost = n->get_cost();
        m_best_config = n->get_assignment();

		return NULL;
	}

    // Keeps track of the highest node to be deleted during cleanup,
    // where .second will be deleted as a child of .first
    std::pair<search_node*, search_node*> highest_delete(NULL, NULL);

    // 'prop' signals whether we are still propagating values in this call
    bool prop = true;
    // 'del' signals whether we are still deleting nodes in this call
    bool del = (upper_limit != n) ? true : false;

    if (m_verbose > 0) {
        std::cout << "[VALUE] begin propagation:" << std::endl;
    }

    // Going all the way to the root, if we have to
    do {

        if (cur->get_type() == MERLIN_NODE_AND) {
            // =================================================================

            // Propagate and update node values
            if (prop) {
                
                // Optimal solution to previously solved and deleted child OR nodes
                double d = cur->get_subsolved();
                
                // Current best solution to yet-unsolved OR child nodes
                std::vector<search_node*>& children = cur->get_children();
                for (size_t i = 0; i < children.size(); ++i) {
                    if (children[i]) {
                        d *= children[i]->get_cost();
                    }
                }

                // store into value (thus includes cost of subSolved)
                cur->set_cost(d);

                if (m_verbose > 0) {
    			    std::cout << "   current AND node updated: " << cur->to_string() << std::endl;
                }

                // Not all OR children solved yet, propagation stops here
                if ( isnan(d) ) { // not all OR children solved yet, propagation stops here
                    prop = false;
                    propagate_tuple(n, cur); // save (partial) opt. subproblem solution at current AND node
                }

            }

            // Clean up fully propagated nodes (i.e. no children or only one (=prev))
            if (del) {
                if (prev->num_children() <= 1) {

                    // prev is OR node, try to cache
                    if (m_caching && prev->is_cachable() && prev->is_optimal()) {
                        try {
                            m_space->write(prev->get_variable(),
                                prev->get_context(), prev->get_cost(),
                                prev->get_assignment() );

                        } catch (...) { /* tried to cache NaN value */
                        }
                    }

                    highest_delete = std::make_pair(cur, prev);
                } else {
                    del = false; // there are unsolved children
                }

            }

            // ===========================================================================
        } else { // cur is OR node
            // ===========================================================================

            if (prop) {
                double d = prev->get_cost() * prev->get_weight(); // getValue includes subSolved

                if (isnan(cur->get_cost()) || d > cur->get_cost()) {
                    cur->set_cost(d); // update max value
                    cur->set_argmax(prev->get_value());
                } else {
                    prop = false; // no more value propagation upwards in this call
                }
            }

            if (m_verbose > 0) {
    			std::cout << "   current OR node updated: " << cur->to_string() << std::endl;
            }

            if (del) {
                if (prev->num_children() <= 1) { // prev has no or one children?
                    highest_delete = std::make_pair(cur, prev);
                } else {
                    del = false;
                }
            }

    		// save opt. tuple, will be needed for caching later
		    if ( prop && cur->is_cachable() && cur->is_optimal() ) {
                //DIAG(myprint("< Cachable OR node found\n"));
                propagate_tuple(n, cur);
            }

            // ===========================================================================
        }

        // Don't delete anything higher than upperLimit
        if (upper_limit == cur) {
            del = false;
        }

        // Move pointers up in search space
        if (prop || del) {
            prev = cur;
            cur = cur->get_parent();
        } else {
            break;
        }

    } while (cur); // until cur==NULL, i.e. 'parent' of root

    if (m_verbose > 0) {
        std::cout << "[VALUE] end propagation." << std::endl;
    }
    
    // propagated up to root node, update tuple as well
    if (prop && !cur) {
        propagate_tuple(n, prev);
        if (solution) {
            double timestamp = timeSystem() - m_start_time;
            m_solutions.push_back(prev->get_assignment());
            update_solution(timestamp,
                    prev->get_cost(),
                    prev->get_assignment(),
                    m_space->get_num_nodes());
        }
    }

    if (highest_delete.first) {
        search_node* parent = highest_delete.first;
        search_node* child = highest_delete.second;
        if (parent->get_type() == MERLIN_NODE_AND) {
            // Store value of OR node to be deleted into AND parent
            parent->add_subsolved(child->get_cost());
        }

        // finally clean up, delete subproblem with unnecessary nodes from memory
        parent->remove_child(child);
        if (m_verbose > 0) {
            std::cout << "[MEMORY] Deleting from memory: " << child->to_string() << std::endl;
        }
        delete child;
    }

    return highest_delete.first;
}

/* collects the joint assignment from 'start' upwards until 'end' and
 * records it into 'end' for later use */
void bound_propagator::propagate_tuple(search_node* start, search_node* end) {

    // Safety checks
	assert(start && end);

    // get the subproblem vars for end node
	int end_var = end->get_variable();
	const std::set<size_t>& end_subprob = m_pseudotree->get_node(end_var)->get_subproblem();

  	// get variable map for end node
	std::vector<int> end_var_map = m_pseudotree->get_node(end_var)->get_subproblem_map();

    // Allocate assignment in end node
	std::vector<int>& assig = end->get_assignment();
	assig.resize(end_subprob.size(), UNKNOWN);

    if (m_verbose > 0) {
        std::cout << "[TUPLE] begin tuple propagation:" << std::endl;
        std::cout << "     end: " << end->to_string() << std::endl;
        std::cout << "   start: " << start->to_string() << std::endl;
        std::cout << " subprob: "; 
        std::copy(end_subprob.begin(), end_subprob.end(), std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "     map: ";
        std::copy(end_var_map.begin(), end_var_map.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
    }

	int curr_var = UNKNOWN, curr_val = UNKNOWN;
	for (search_node* curr = start; curr != end; curr = curr->get_parent()) {
		curr_var = curr->get_variable();

        if (m_verbose) {
            std::cout << "  > at: " << curr->to_string() << std::endl;
        }

		if (curr->get_type() == MERLIN_NODE_AND) {
			curr_val = curr->get_value();
			if (curr_val != UNKNOWN) {
			    assig.at(end_var_map.at(curr_var)) = curr_val;
            }

            if (m_verbose > 0) {
                std::cout << "assg: ";
                std::copy(assig.begin(), assig.end(), std::ostream_iterator<int>(std::cout, " "));
                std::cout << std::endl;
            }
		}

		if (curr->get_assignment().size() > 0) {
			// check previously saved partial assignment
			const std::set<size_t>& curr_subprob = m_pseudotree->get_node(curr->get_variable())->get_subproblem();
			std::set<size_t>::const_iterator itVar = curr_subprob.begin();
			std::vector<int>::const_iterator itVal = curr->get_assignment().begin();

			for(; itVar!= curr_subprob.end(); ++itVar, ++itVal ) {
                size_t var = *itVar;
                int val = *itVal;
                if (*itVal != UNKNOWN) {
                    assig[end_var_map[var]] = val;
                }
			}

			// clear optimal assignment of AND node, since now propagated upwards
			if (curr->get_type() == MERLIN_NODE_AND) {
			    curr->clear_assignment();// TODO correct ?
            }

            if (m_verbose > 0) {
                std::cout << "assg: ";
                std::copy(assig.begin(), assig.end(), std::ostream_iterator<int>(std::cout, " "));
                std::cout << std::endl;
            }
		}
	} // end for

    if (m_verbose > 0) {
        std::cout << "[TUPLE] end tuple propagation." << std::endl;
    }
}

void bound_propagator::update_solution(double timestamp, double cost, 
    std::vector<int>& sol, std::pair<size_t, size_t> num_nodes) {

    if ( (isnan(m_best_cost) || cost > m_best_cost) ) {
        m_best_cost = cost;
        m_best_config = sol;
        std::cout << "[" << std::setw(9) << timestamp << "] u "
            << std::setw(12) << num_nodes.first << " "
            << std::setw(12) << num_nodes.second << " "
            << std::setw(12) << m_best_cost << " (" << std::log10(m_best_cost) << ")"
            << std::endl;
    } else {
        return;
    }
}

} // end namespace


 