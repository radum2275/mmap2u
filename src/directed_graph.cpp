/*
 * directed_graph.cpp
 *
 *  Created on: Oct 06, 2020
 *      Author: radu
 *
 * Copyright (c) 2020, International Business Machines Corporation
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

#include "directed_graph.h"
#include "graph.h"
#include "utils.h"

namespace merlin {

const directed_edge directed_edge::NO_EDGE(-1, -1, 0);

std::vector<size_t> directed_graph::topological_sort() {
	std::vector<size_t> indegree(num_nodes(), 0);
	for (size_t i = 0; i < num_nodes(); ++i) {
		const my_set<directed_edge>& n = neighbors(i);
		my_set<directed_edge>::const_iterator it = n.begin();
		for (; it != n.end(); ++it) {
			assert(it->first == i);
			indegree[it->second]++;
		}
	}

	std::vector<size_t> result;
	std::queue<size_t> Q;
	for (size_t i = 0; i < indegree.size(); ++i) {
		if (indegree[i] == 0) {
			Q.push(i);
		}
	}

	while (Q.empty() == false) {
		size_t v = Q.front();
		Q.pop();

		assert(indegree[v] == 0);
		result.push_back(v);

		const my_set<directed_edge>& n = neighbors(v);
		my_set<directed_edge>::const_iterator it = n.begin();
		for (; it != n.end(); ++it) {
			size_t u = it->second;
			indegree[u]--;
			if (indegree[u] == 0) {
				Q.push(u);
			}
		}
	}

	return result;
}

// Check if it is singly connected
bool directed_graph::is_cyclic() {
	size_t n = num_nodes();
	graph g(n); // create the undirected graph
	for (size_t v = 0; v < n; ++v) {
		for (my_set<directed_edge>::const_iterator ei = m_adj[v].begin(); 
			ei != m_adj[v].end(); ++ei) {
			g.add_edge(ei->first, ei->second);
		}
	}

	bool result = g.is_cyclic();
	return result;
}

// Find a loop cutset (greedy)
my_set<directed_edge> directed_graph::find_loop_cutset() {
	my_set<directed_edge> cutset;
	directed_graph g(*this);
	bool found = false;

	while (!found) {
		// Get the nodes and their out-degrees
		std::vector<std::pair<size_t, size_t> > nodes;
		for (size_t v = 0; v < g.num_nodes(); ++v) {
			nodes.push_back(std::make_pair(v, g.neighbors(v).size()));
		}

		// Sort the nodes in decreasing order of their out-degrees
		std::sort(nodes.begin(), nodes.end(),
			[](const std::pair<size_t, size_t>& left, std::pair<size_t,size_t>& right) {
				return left.second > right.second;
			}
		);

		// Select the node with largest out-degree (bracking ties randomly)
		std::vector<size_t> candidates;
		size_t max_degree = nodes[0].second;
		if (max_degree == 1) {
			found = true;
		} else {
			candidates.push_back(nodes[0].first);
			size_t i = 1;
			while (nodes[i].second == max_degree) {
				candidates.push_back(nodes[i].first);
				i++;
			}
			std::random_shuffle(candidates.begin(), candidates.end());
			size_t u = candidates[0];

			// Select randomly an outgoing edge
			int n = (int)g.neighbors(u).size();
			size_t j = (size_t) randi(n);
			const directed_edge& e = g.neighbors(u).at(j);
			cutset |= e;
			g.remove_edge(e.first, e.second);

			if (g.is_cyclic() == false) {
				found = true;
			}
		}
	}

	return cutset;
}


} // end namespace
