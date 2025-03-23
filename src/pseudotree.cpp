/*
 * pseudotree.h
 *
 *  Created on: 17 Mar 2025
 *      Author: radu
 *
 * Copyright (c) 2025, International Business Machines Corporation
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

/// \file pseudotree.h
/// \brief Pseudo tree
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#include "pseudotree.h"

namespace merlin {


// Build the pseudo tree
void pseudotree::build(graph& g, std::vector<size_t>& order, bool is_chain) {
    
    // Save the elimination order
    m_order = order;

    // Triangulate the graph
    graph gg(g);
    gg.triangulate(order);

    // std::cout << "[Triangulated moral graph]" << std::endl;
    // const std::vector<edge_id>& edges = gg.edges();
    // for (size_t i = 0; i < edges.size(); ++i) {
    //     std::cout << edges[i] << std::endl;
    // }

    std::list<pseudotree_node*> roots;
    std::cout << "[Building the pseudo tree...]" << std::endl;
    if (is_chain) { // Build a chain pseudo tree (OR search)

        // Search is done in reversed elimination order
        std::vector<size_t> search_order(order.rbegin(), order.rend());
        std::vector<size_t>::const_iterator li = search_order.begin();
        pseudotree_node* r = m_nodes[*li].get();
        r->set_parent(NULL);
        pseudotree_node* prev = r;
        li++;
        while ( li != search_order.end() ){
            size_t var = (*li);
            pseudotree_node* n = m_nodes[var].get();
            n->set_parent(prev);
            prev->add_child(n);
            prev = n;
            li++;
        }

        // Save the current root
        roots.push_back(r);
    } else { // Build the pseudo tree (AND/OR search)

        for (std::vector<size_t>::iterator it = order.begin(); it != order.end(); ++it) {
            size_t v = (*it);
            std::cout << "...processing var:";
            std::set<size_t> N = gg.get_neighbors(v); // set of neighboring nodes
            std::copy(N.begin(), N.end(), std::ostream_iterator<size_t>(std::cout, " "));
            std::cout << std::endl;
            m_width = std::max(m_width, N.size());
            insert_node(v, N, roots);
            gg.add_clique(N);
            gg.remove_node(v, true);
        }
    }

    std::cout << "[Pseudo tree built.]" << std::endl;
	std::cout << "[Connecting disconnected components: " << roots.size() << "]" << std::endl;

    // Add artificial root node to connect disconnected components
	int dummy = order.size();
    m_nodes.push_back(std::make_unique<pseudotree_node>(dummy));
    pseudotree_node* p = m_nodes.back().get();
    std::list<pseudotree_node*>::iterator it = roots.begin();
	for (; it != roots.end(); ++it) {
		pseudotree_node* n = (*it);
        p->add_child(n);
		n->set_parent(p);
	}

    // Set the root of the pseudo tree
	m_root = p;

    // Update the subproblem variables (recursively)
    m_root->update_subproblem(m_nodes.size()); // +1 dummy variable

    // Update contexts if chain pseudo tree
    if (is_chain) {
        update_contexts(gg);
    }

	// Update the elimination order
	m_order.push_back(dummy); // add dummy variable as last node in ordering
}

// Reset the potentials associated with each variable
void pseudotree::reset_potentials(std::vector<interval>& factors) {
    
    // Reset the previous mapping
    m_potentials.clear();
    m_potentials.resize(m_nodes.size());

    std::vector<interval>::iterator vi = factors.begin();
	for (; vi != factors.end(); ++vi) {
        interval& f = (*vi);
        const variable_set& scope = f.vars();
		if (scope.size() == 0) {
            size_t v = m_order.back();
			m_potentials[v].push_back(f.to_potential(false));
			continue;
		}

		std::vector<size_t>::iterator it = m_order.begin();
		for (;; ++it) {
            size_t v = (*it);
			if (scope.has_variable(v)) {
				m_potentials[v].push_back(f.to_potential(false));
				break;
			}
		}
	}

}

// Update the contexts of the nodes
void pseudotree::update_contexts(graph& g) {

	for (size_t i = 0; i < m_nodes.size(); ++i) {

		pseudotree_node *n = m_nodes[i].get();
		if (n->get_parent() == NULL) {
			continue; // skip the root
		}

		std::set<size_t> ancestors; // this is the context
		const std::set<size_t>& descendants = n->get_subproblem();

		// Find the ancestors of the node that are connected to it or to its 
        // descendants in the triangulated moral graph (i.e., induced graph)
		pseudotree_node *p = n->get_parent();
		while (p != NULL) {

			bool found = false;
			for (std::set<size_t>::const_iterator si = descendants.begin();
					si != descendants.end(); ++si) {
				
                size_t d = (*si); // includes current node as well
                size_t v = p->get_variable();

				if (g.edge(v, d) != edge_id::NO_EDGE) {
					found = true;
					break;
				}
			}

			// Found an ancestor connected to a descendant in the induced graph
			if (found) {
				ancestors.insert(p->get_variable());
			}

			p = p->get_parent();
		}

		// Set the context
		n->set_context(ancestors);
	}

}

// Insert a new node in the pseudo tree
void pseudotree::insert_node(size_t v, std::set<size_t>& neighbors, 
    std::list<pseudotree_node*>& roots) {

    // Safety checks
    assert(v >= 0 && v < m_nodes.size());

    // Create new node in pseudo tree
    pseudotree_node* p = m_nodes[v].get();
    p->set_context(neighbors);

	// Insert the new pseudo tree node
	std::list<pseudotree_node*>::iterator it = roots.begin();
	while (it != roots.end()) {
        pseudotree_node* n = (*it);
		const std::set<size_t>& context = n->get_context();
		if (context.find(v) != context.end()) {
			p->add_child(n); // add child to current node
			n->set_parent(p); // set parent of previous node
			it = roots.erase(it); // remove previous node from roots list
		} else {
			++it;
		}
	}

	roots.push_back(p); // add current node to roots list
}

// Dump the pseudo tree to an output stream
void pseudotree::dump(std::ostream& os) {
    os << "Pseudo tree nodes: " << m_nodes.size() << std::endl;
    os << "Pseudo tree root : " << m_root->get_variable() << std::endl;
    for (size_t v = 0; v < m_nodes.size(); ++v) {
        pseudotree_node* n = m_nodes[v].get();
        pseudotree_node* p = n->get_parent();
        os << "[" << n->get_variable() << "] - ";
        if (p == NULL) {
            os << "():"; 
        } else {
            os << "(" << p->get_variable() << "):";
        }

        const std::vector<pseudotree_node*>& children = n->get_children();
        for (size_t i = 0; i < children.size(); ++i) {
            os << " " << children[i]->get_variable();
        }

        os << std::endl;
    }

    os << "Subproblems:" << std::endl;
    for (size_t v = 0; v < m_nodes.size(); ++v) {
        pseudotree_node* n = m_nodes[v].get();
        os << "[" << n->get_variable() << "]:";
        const std::set<size_t>& subproblem = n->get_subproblem();
        std::copy(subproblem.begin(), subproblem.end(), std::ostream_iterator<size_t>(os, " "));
        os << std::endl;
    }

    os << "Contexts:" << std::endl;
    for (size_t v = 0; v < m_nodes.size(); ++v) {
        pseudotree_node* n = m_nodes[v].get();
        os << "[" << n->get_variable() << "]:";
        const std::set<size_t>& context = n->get_context();
        std::copy(context.begin(), context.end(), std::ostream_iterator<size_t>(os, " "));
        os << std::endl;
    }

}

// Dump the pseudo tree to a Graphviz dot file
void pseudotree::dump_for_dot(std::string filename) {
    std::ofstream outfile(filename.c_str());
    if (outfile.fail()) {
        std::cerr << "Cannot open Graphviz file.";
        exit(EXIT_FAILURE);
    }

    outfile << "digraph g {\n";
    outfile << "node [shape = record];\n";
    outfile << "size = \"10, 7.5\";\n";
    outfile << "rotate = \"90\";\n";
    outfile << "ratio = \"fill\";\n";
    dump_nodes_for_dot(outfile);
    dump_edges_for_dot(outfile);
    outfile << "}\n";
}

// Dump the nodes in the pseudo tree.
void pseudotree::dump_nodes_for_dot(std::ofstream& outfile ) {

	std::queue<pseudotree_node*> q;

	q.push(m_root);
	while (!q.empty()) {
		pseudotree_node* node = q.front();
		q.pop();

    	outfile << "node" << node->get_variable()
			<< "[ shape=ellipse, color=gold, label = \"" << node->get_variable() << "\"];\n";

		const std::vector<pseudotree_node*>& ch = node->get_children();
		for (size_t i = 0; i < ch.size(); ++i) {
			pseudotree_node* child = ch[i];
			q.push(child);
		}
	}
}

// Dump the edges in the pseudo tree.
void pseudotree::dump_edges_for_dot(std::ofstream& outfile) {
	
    std::queue<pseudotree_node*> q;

	q.push(m_root);
	while (!q.empty()) {
		pseudotree_node* node = q.front();
		q.pop();

		const std::vector<pseudotree_node*>& ch = node->get_children();
		for (size_t i = 0; i < ch.size(); ++i) {

			pseudotree_node* child = ch[i];

			outfile << "node" << node->get_variable()
				<< " -> node" << child->get_variable() << ";\n";

			q.push(child);
		}
	}
}

// Update the subproblem rooted at the node
const std::set<size_t>& pseudotree_node::update_subproblem(size_t num_vars) {

	// Clear current subproblem
	m_subproblem.clear();
	// Add self (i.e., node variable)
	m_subproblem.insert(m_variable);

	// Iterate over children and collect their subproblem variables
	for (std::vector<pseudotree_node*>::iterator it = m_children.begin();
			it != m_children.end(); ++it) {
        
        pseudotree_node* ch = (*it);
		const std::set<size_t>& child_vars = ch->update_subproblem(num_vars);
        std::copy(child_vars.begin(), child_vars.end(), 
            std::inserter(m_subproblem, m_subproblem.end()));
	}

    m_subproblem_map.clear();
	m_subproblem_map.resize(num_vars, UNKNOWN);
	size_t i = 0;
	for (std::set<size_t>::const_iterator it = m_subproblem.begin();
			it != m_subproblem.end(); ++it, ++i) {
		m_subproblem_map[*it] = i;
	}

	// Return a const reference
	return m_subproblem;
}

} // end namespace

