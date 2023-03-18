/*
 * graph.cpp
 *
 *  Created on: 15 May 2015
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

#include "graph.h"
#include "utils.h"

namespace merlin {

const edge_id edge_id::NO_EDGE(-1, -1, 0, 0);

/// Triangulate the graph along an ordering
void graph::triangulate(const std::vector<size_t>& ordering) {

	// set the position map
	size_t n = ordering.size();
	std::vector<size_t> position(n);
	for (size_t i = 0; i < n; ++i) {
		position[ordering[i]] = i;
	}

	// triangulate the graph
	for (size_t i = 0; i < n; ++i) {
		size_t var = ordering[i];
		size_t pos = position[var];

		const my_set<edge_id>& ns = neighbors(var);
		std::vector<size_t> tmp;
		for (my_set<edge_id>::const_iterator si = ns.begin(); si != ns.end(); ++si) {
			index ni = (*si).second;
			if (position[ni] > pos) {
				tmp.push_back(ni);
			}
		}

		// connect the earlier neighbors
		int m = (int)tmp.size();
		for (int j = 0; j < m-1; ++j) {
			for (int k = j+1; k < m; ++k) {
				add_edge(tmp[j], tmp[k]);
			}
		}
	}
}


/// Create a graph from a set of factor scopes
void graph::init(const my_vector<variable_set>& fin) {
	for (my_vector<variable_set>::const_iterator i = fin.begin();
			i != fin.end(); ++i) {
		const variable_set& vs = (*i);
		for (size_t ii = 0; ii < vs.size()-1; ++ii) {
			for (size_t jj = ii+1; jj < vs.size(); ++jj) {
				add_edge(vs[ii].label(), vs[jj].label());
			}
		}
	}
}

void graph::init(const std::vector<std::vector<size_t> >& scopes) {
	for (size_t i = 0; i < scopes.size(); ++i) {
		const std::vector<size_t>& vs = scopes[i];
		for (size_t ii = 0; ii < vs.size() - 1; ++ii) {
			for (size_t jj = ii+1; jj < vs.size(); ++jj) {
				add_edge(vs[ii], vs[jj]);
			}
		}
	}
}

// Check if the graph is disconnected. Run a DFS traversal and check if all 
// nodes were visited. If not, then the graph is disconnected.
bool graph::is_disconnected() {
	size_t n = num_nodes();
	std::vector<bool> visited(n, false);

	std::stack<size_t> dfs;
	dfs.push(0);
	while (!dfs.empty()) {
		size_t vid = dfs.top();
		dfs.pop();
		visited[vid] = true;
		
		const my_set<edge_id>& ns = neighbors(vid);
		for (my_set<edge_id>::const_iterator si = ns.begin();
				si != ns.end(); ++si) {
			size_t ni = (*si).second;
			if (!visited[ni]) {
				dfs.push(ni);
			}
		}
	}
	
	for (size_t i = 0; i < visited.size(); ++i) {
		if (!visited[i]) {
			return true;
		}
	}

	return false;
}

/// Retrieve the cliques of a triangulated graph
std::vector<std::set<size_t> > graph::maximal_cliques(const std::vector<size_t>& ordering) {

	// Create the position map
	size_t n = ordering.size();
	std::vector<size_t> position(n);
	for (size_t i = 0; i < n; ++i) {
		position[ordering[i]] = i;
	}

	// Get all cliques
	std::vector<std::set<size_t> > cliques;
	for (size_t i = 0; i < ordering.size(); ++i) {
		size_t var = ordering[i];
		size_t pos = position[var];
		const my_set<edge_id>& ns = neighbors(var);

		// form the clique
		std::set<size_t> cli;
		cli.insert(var);
		for (my_set<edge_id>::const_iterator si = ns.begin();
				si != ns.end(); ++si) {
			size_t ni = (*si).second;
			if (position[ni] > pos) {
				cli.insert(ni);
			}
		}

		cliques.push_back(cli);
	}

	// Keep only maximal cliques
	// remove dominated cliques
	// : A dominates B ( A > B ) iff scope(A) includes scope(B)
	std::vector<std::set<size_t> > clusters;
	for (std::vector<std::set<size_t> >::iterator it1 = cliques.begin();
			it1 != cliques.end(); ++it1) {

		std::set<size_t>& A = *it1;
		bool found = false;
		for (std::vector<std::set<size_t> >::iterator it2 = clusters.begin();
				it2 != clusters.end(); ++it2) {

			std::set<size_t>& B = *it2;
			if ( dominates(B, A) ) {
				found = true;
				break;
			}
		}

		if (!found) {
			for (std::vector<std::set<size_t> >::iterator it3 = clusters.begin();
					it3 != clusters.end();) {

				std::set<size_t>& B = (*it3);
				if ( dominates(A, B) ) {
					clusters.erase(it3);
				} else {
					++it3;
				}
			}

			clusters.push_back(A);
		}
	}

	return clusters;
}

bool graph::cyclic_util(size_t v, std::vector<bool>& visited, int parent) { 
      
    // Mark the current node as visited 
    visited[v] = true; 
  
    // Recur for all the vertices adjacent to this vertex 
	const my_set<edge_id>& neigh = neighbors(v);
	my_set<edge_id>::const_iterator i;
    for (i = neigh.begin(); i != neigh.end(); ++i) { 
        assert(i->first == v);
		size_t u = i->second;
        // If an adjacent is not visited, then recur for that adjacent 
        if (!visited[u]) {
			if (cyclic_util(u, visited, v)) {
				return true; 
			}
        } else if ((int)u != parent) {
			// If an adjacent is visited and not parent of current vertex, then there is a cycle. 
			return true; 
		}
    }

    return false; 
} 
  
// Returns true if the graph contains a cycle, else false. 
bool graph::is_cyclic() { 
      
    // Mark all the vertices as not visited and not part of recursion stack 
    std::vector<bool> visited(num_nodes(), false); 
  
    // Call the recursive helper function to detect cycle in different DFS trees 
    for (size_t u = 0; u < num_nodes(); ++u) { 
        // Don't recur for u if it is already visited 
        if (!visited[u])
			if (cyclic_util(u, visited, -1)) 
            	return true; 
    }

    return false; 
} 

} // end namespace
