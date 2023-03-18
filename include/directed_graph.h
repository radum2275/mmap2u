/*
 * directed_graph.h
 *
 *  Created on: Oct 6, 2020
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

/// \file directed_graph.h
/// \brief An directed graph structure
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_DIRECTED_GRAPH_H_
#define IBM_LOOPY_DIRECTED_GRAPH_H_

#include "base.h"
#include "set.h"

namespace merlin {

typedef std::pair<size_t, size_t> edge_t; ///< Basic Edge type.

///
/// \brief Edge in the graph.
///
/// EdgeID represents the edge (i->j) with index *idx*. Blank (or missing) edges 
/// can sometimes appear, and they are represented by EdgeID::NO_EDGE 
/// (for example, if edge (i,j) is  requested and does not exist).
///
struct directed_edge {
	typedef size_t index; 	///< Basic indexing type for the edge.
	index first;			///< Head of the edge.
	index second;			///< Tail of the edge.
	index idx;				///< Index of the edge.

	///
	/// \brief Default constructor.
	///
	directed_edge() : first(), second(), idx() {};

	///
	/// \brief Constructor.
	///
	/// Create the edge (i->j) in the graph.
	/// \param i 	Head of the edge
	/// \param j 	Tail of the edge
	/// \param eij 	Index of the edge
	/// \param eji 	Index of the reversed edge
	directed_edge(index i, index j, index eij) 
		: first(i), second(j), idx(eij) {};

	///
	/// \brief Indexing operator.
	///
	operator edge_t() const { 
		return edge_t(first,second); 
	}

	static const directed_edge NO_EDGE; ///< A static constant for representing missing edges.

	///
	/// \brief Comparison (equal) operator.
	///
	bool operator==(const directed_edge& e) const { 
		return e.first==first && e.second==second; 
	};
	
	///
	/// \brief Comparison (not equal) operator.
	///
	bool operator!=(const directed_edge& e) const { 
		return !(*this==e); 
	};
	
	///
	/// \brief Comparison (less than) operator.
	///
	bool operator< (const directed_edge& e) const { 
		return (first < e.first || (first==e.first && second<e.second)); 
	};
};

///
/// \brief The graph structure used by Merlin.
///
/// The graph contains N *nodes* (indexed 0..N-1) and E directed *edges*, with
/// indeces (0..E-1).
///
class directed_graph {
public:
	typedef directed_edge::index index;			///< Basic indexing type for the graph.
	std::vector<my_set<directed_edge> > m_adj;		///< Look up EdgeID info by adj[i][jj] (jj = position of j).
	std::vector<directed_edge> m_edges;			///< Look up EdgeID info by edges[eij] (edge index).

protected:
	std::stack<index> m_vvacant;		///< List of available vertex ids.
	std::stack<index> m_evacant;		///< List of available edge ids.

public:
	///
	/// \brief Constructs a graph with *n* vertices and 0 edges.
	/// \param n 	The number of nodes
	///
	explicit directed_graph(size_t n=0) : m_adj(), m_edges(), m_vvacant(), m_evacant() {
		m_adj.resize(n);
	};

	///
	/// \brief Copy constructor
	///
	directed_graph(const directed_graph& g) : m_adj(g.m_adj), m_edges(g.m_edges),
		m_vvacant(g.m_vvacant), m_evacant(g.m_evacant) {};

	///
	/// \brief Destroys the graph.
	///
	~directed_graph() {
	};

	///
	/// \brief Add a node to the graph.
	/// \return the index of the newly added node.
	///
	index add_node() {			// verify that the adjacency table is large enough
		index use;				// and get an index off the stack if available
		if (m_vvacant.empty()) {
			use = num_nodes();
			m_adj.resize(num_nodes()+1);
		} else {
			use = m_vvacant.top(); m_vvacant.pop();
		}
		return use;
	};

	///
	/// \brief Add the nodes to the graph.
	///
	void add_nodes(size_t n) {
		m_adj.resize(n);
	}

	///
	/// \brief Remove a node from the graph.
	/// \param i 	Index of the node to be removed
	///
	void remove_node(index i) {		// pop_back, or add to vVacant
		m_vvacant.push(i); 	// can we keep track of "real" nNodes and iterate through them?
	};

	///
	/// \brief Return the number of nodes.
	///
	size_t num_nodes() const {
		return m_adj.size();
	};

	///
	/// \brief Return the number of edges.
	///
	size_t num_edges() const {
		return m_edges.size();
	};

	///
	/// \brief Add an edge to the graph. 
	/// Create the pairs (i,j) in the adjacency list.
	/// \param i 	The head of the edge
	/// \param j 	The tail of the edge
	/// \return the id (index) of the newly added edge.
	///
	const directed_edge& add_edge(index i, index j) { // add edge (i,j) to adj
		//std::cout<<"Add edge "<<i<<","<<j<<"\n";
		if (edge(i,j) != directed_edge::NO_EDGE)
			return edge(i,j);	// if exists already do nothing

		size_t eij, emax = num_edges();	// otherwise get two edge indices
		if (m_evacant.empty()) {
			eij = emax++;
			m_edges.resize(emax);
		} else {
			eij = m_evacant.top();
			m_evacant.pop();
		}

		m_edges[eij] = directed_edge(i, j, eij);
		m_adj[i] |= m_edges[eij];

		return m_edges[eij];
	};

	///
	/// \brief Remove an edge from the graph.
	/// \param i 	The head of the edge to be removed
	/// \param j 	The tail of the edge to be removed
	///
	void remove_edge(index i, index j) {	// remove edge (i,j) from adj; add to m_evacant
		if (edge(i,j) != directed_edge::NO_EDGE) {
			directed_edge e = edge(i,j);
			index eij = e.idx;
			m_evacant.push(eij); // add eij to edge stack
			m_adj[i] /= m_edges[eij];
			m_edges[eij]=directed_edge::NO_EDGE;
		}
	};

	///
	/// \brief Clear the graph by removing all nodes and edges.
	///
	void clear() {
		clear_edges();
		while (!m_vvacant.empty()) m_vvacant.pop();
		m_adj.resize(0);
	};

	///
	/// \brief Remove all edges.
	///
	void clear_edges() {
		m_edges.clear();
		size_t n = num_nodes(); m_adj.clear(); m_adj.resize(n);
		while (!m_evacant.empty()) m_evacant.pop();
	};

	///
	/// \brief Return the edge (if any) between two nodes.
	/// \param i 	The head of the edge
	/// \param j 	The tail of the edge
	/// \return the index of the edge if the edge is present or 
	/// 	the static constant NO_EDGE if the edge is missing.
	///
	const directed_edge& edge(index i, index j) const {
		if (i >= m_adj.size()) return directed_edge::NO_EDGE;
		my_set<directed_edge>::const_iterator it = m_adj[i].find( directed_edge(i, j, 0) );
		if (it==m_adj[i].end()) return directed_edge::NO_EDGE;
		else return *it;
	};

	///
	/// \brief Return the neighbors of a node (ie, adjacency list).
	/// \param i 	The index of the node
	/// \return the set of edge ids that are adjacent to the node.
	///
	const my_set<directed_edge>& neighbors(index i) const {
		return m_adj[i];
	};

	///
	/// \brief Return an edge (by its index).
	/// \param eij 	The index of the edge
	/// \return the edge id corresponding to the index.
	const directed_edge& edge(index eij) const {
		return m_edges[eij];
	};

	///
	/// \brief Return the edges of the graph.
	/// \return the list of all edge ids in the graph.
	///
	const std::vector<directed_edge>& edges() const {
		return m_edges;
	};

	///
	/// \brief Return a topological order of the vertices
	///
	std::vector<size_t> topological_sort();

	///
	/// \brief Check if the graph is singly connected or not
	///
	bool is_cyclic();

	///
	/// \brief Find a loop cutset (greedy)
	///
	my_set<directed_edge> find_loop_cutset();
};

///
/// \brief Output operator.
///
/// Write out the (formatted) content of an edge.
/// \param out 	The output stream to be written
/// \param e 	The edge to be written out
/// \return a reference to the modified output stream containing the edge
///		contents.
///
inline std::ostream& operator<<( std::ostream& out, const directed_edge& e) {
	return out << "(" << (int)e.first << "," << (int)e.second <<")";
}

} // namespace

#endif /* IBM_LOOPY_DIRECTED_GRAPH_H_ */
