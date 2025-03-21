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

#ifndef IBM_LOOPY_PSEUDOTREE_H_
#define IBM_LOOPY_PSEUDOTREE_H_

#include "base.h"
#include "graph.h"
#include "interval.h"
#include "potential.h"

namespace merlin {

class pseudotree_node;

///
/// @brief Pseudo tree structure
///
class pseudotree {
public:
    typedef std::unique_ptr<pseudotree_node> pseudotree_node_ptr;

public:
    /// @brief Constructor.
    pseudotree() : m_root(NULL), m_height(0), m_width(0) {};

    /// @brief Destructor.
    ~pseudotree() {
        m_nodes.clear();
    }

    /// @brief Initialize the pseudo tree.
    /// @param n is the number of variables
    void init(size_t n) {
        m_nodes.clear();
        for (size_t v = 0; v < n; ++v) {
            m_nodes.push_back(std::make_unique<pseudotree_node>(v));
        }    
    }

    /// @brief Get the root of the pseudo tree.
    /// @return A pointer to the root node.
    inline pseudotree_node* get_root(){
        return m_root;
    }

    /// @brief Build the pseudo tree from a directed graph.
    /// @param g is a directed graph
    /// @param order is the elimination order
    /// @param is_chain is a flag indicating a chain pseudo tree (OR search)
    void build(graph& g, std::vector<size_t>& order, bool is_chain = false);

    /// @brief Dump the pseudo tree in a Graphviz dot file.
    /// @param filename is the name of the output dot file
    void dump_for_dot(std::string filename);

    /// @brief Dump the pseudo tree to an output stream.
    /// @param os is the output stream
    void dump(std::ostream& os);

    /// @brief Get the pseudo tree node corresponding to a variable index.
    /// @param v is the variable index
    /// @return The pseudo tree node corresponding to the variable.
    pseudotree_node* get_node(size_t v) {
        assert(v >= 0 && v < m_nodes.size());
        return m_nodes[v].get();
    }

    /// @brief Get the induced width of the pseudo tree.
    /// @return The induced width.
    inline size_t get_width() {
        return m_width;
    }

    /// @brief Get the height of the pseudo tree.
    /// @return The height.
    inline size_t get_height() {
        return m_height;
    }
    
    /// @brief Get the original potentials associated with a variable
    /// @param var is the index of the variable
    /// @return a list of potentials.
    std::list<potential>& get_potentials(size_t var) {
        return m_potentials[var];
    }

   /// @brief Reset the mapping of potentials to variabeles.
   void reset_potentials(std::vector<interval>& factors);

protected:

    /// @brief Add a new node to the pseudo tree.
    /// @param v is the variable index associated with the new node
    /// @param neighbors is the set of neighbors of v in the graph
    /// @param roots is the temporary list of roots of the pseudo tree
    void insert_node(size_t v, std::set<size_t>& neighbors, std::list<pseudotree_node*>& roots);

    /// @brief Dump the pseudo tree nodes into a dot file.
    /// @param os is the output file stream
    void dump_nodes_for_dot(std::ofstream& os);

    /// @brief Dump the pseudo tree edges into a dot file.
    /// @param os is the output file stream
    void dump_edges_for_dot(std::ofstream& os);

    /// @brief Update the contexts of the nodes
    /// @param g is the triangulated moral graph
    void update_contexts(graph& g);

protected:
    // Members

    std::vector<pseudotree_node_ptr> m_nodes;       ///< Nodes of the pseudo tree
    size_t m_height;                                ///< Height of the pseudo tree
    size_t m_width;                                 ///< Width of the pseudo tree
    std::vector<size_t> m_order;                    ///< Elimination order
    pseudotree_node* m_root;                        ///< The root node of the pseudo tree
    std::vector<std::list<potential>> m_potentials; ///< The original potentials

};

/// @brief Pseudo tree node structure
class pseudotree_node {
public:

    /// @brief Constructor.
    /// @param v is the variable index of the node
    pseudotree_node(size_t v)
        : m_variable(v), m_parent(NULL) {};

    /// @brief Constructor.
    /// @param v is the variable index of the node
    /// @param neighbors is the set of neighboring variables in the graph (OR context)
    pseudotree_node(size_t v, std::set<size_t>& neighbors)
        : m_variable(v), m_context(neighbors), m_parent(NULL) {};

    /// @brief Destructor.
    ~pseudotree_node() {};


    /// @brief Get the node variable.
    /// @return The variable index.
    inline size_t get_variable() {
        return m_variable;
    }

    /// @brief Set the parent of the current node.
    /// @param p is the parent node inthe pseudo tree
    inline void set_parent(pseudotree_node* p) {
        m_parent = p;
    }

    inline pseudotree_node* get_parent() {
        return m_parent;
    }

    /// @brief Add a child to the current node.
    /// @param c is the child node in the pseudo tree
    inline void add_child(pseudotree_node* c) {
        m_children.push_back(c);
    }

    /// @brief Get the children of a node.
    /// @return A const reference to the vector containing the children.
    inline const std::vector<pseudotree_node*>& get_children() const {
        return m_children;
    }

    /// @brief Get the context of a node in the pseudo tree.
    /// @return A const reference to the set representing the context.
    const std::set<size_t>& get_context() const {
        return m_context;
    }

    /// @brief Update the subproblem rooted at the node (includes the node)
    /// @return The subproblem variables.
    const std::set<size_t>& update_subproblem();

    /// @brief Get the subproblem variables.
    /// @return A const reference to the subproblem variables.
    inline const std::set<size_t>& get_subproblem(){
        return m_subproblem;
    }

    /// @brief Set the context of a node in the pseudo tree (does not include variable).
    /// @param context is the set of variables representing the context (OR contexts)
    void set_context(std::set<size_t>& context) {
        m_context = context;
    }

    /// @brief Dump the node context.
    /// @param os is the output stream
    inline void dump_contexts(std::ostream& os) {
        os << m_variable << ": [ ";
        std::copy(m_context.begin(), m_context.end(),
                std::ostream_iterator<size_t>(os, " "));
        os << "]" << std::endl;
    }

protected:
    // Members

    size_t m_variable;                          ///< The variable index
    std::set<size_t> m_context;                 ///< The context (OR)
    pseudotree_node* m_parent;                  ///< The parent node in the pseudo tree
    std::vector<pseudotree_node*> m_children;   ///< The children nodes in the pseudo tree
    std::set<size_t> m_subproblem;              ///< The subproblem rooted at the node (includes the node)

};

} // end namespace

#endif
