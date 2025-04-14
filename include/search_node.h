/*
 * search_node.h
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

/// \file search_node.h
/// \brief Search node
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_SEARCH_NODE_H_
#define IBM_LOOPY_SEARCH_NODE_H_

#include "base.h"

namespace merlin {

///
/// @brief A search node structure
///
class search_node {
    protected:
        size_t m_variable;						// variable
        int m_value;							// value
        size_t m_type;							// node type (AND, OR)
        double m_heur;							// as in h-value
        double m_cost;							// as in g-value
        double m_weight;                        // OR-AND arc weight
        size_t m_depth;                         // depth
        bool m_leaf;                            // leaf node
        bool m_cachable;                        // cachable node
        bool m_pruned;                          // pruned node
        bool m_optimal;                         // optimal node (true)
        double m_subsolved;                     // cost of solved subproblems
        int m_argmax;                           // best domain value for OR nodes

        search_node* m_parent;                  // parent
        std::vector<search_node*> m_children;   // children
        std::vector<double> m_cache;            // heuristic cache
        std::vector<int> m_assignment;          // optimal assignment below node
        std::string m_context;                  // cache context

    public:
        search_node(size_t var, int val, size_t type) :
            m_variable(var), m_value(val), m_type(type), m_heur(1.0), 
            m_cost(NAN), m_weight(1.0), m_parent(NULL), m_leaf(false),
            m_subsolved(1.0), m_cachable(false), m_pruned(false), 
            m_optimal(true), m_argmax(UNKNOWN) {};
        
        ~search_node() {};

        inline size_t get_variable() {
            return m_variable;
        }
        inline int get_value() {
            return m_value;
        }
        inline size_t get_type() {
            return m_type;
        }
        inline double get_heur() const {
            return m_heur;
        }
        inline void set_heur(double h) {
            m_heur = h;
        }
        inline double get_cost() const {
            return m_cost;
        }
        inline double get_weight() {
            return m_weight;
        }
        inline size_t get_depth() {
            return m_depth;
        }
        inline double get_subsolved() {
            return m_subsolved;
        }
        inline void set_subsolved(double d) {
            m_subsolved = d;
        }
        inline void add_subsolved(double d) {
            m_subsolved *= d;
        }
        inline void set_weight(double w) {
            m_weight = w;
        }
        inline void set_depth(size_t d) {
            m_depth = d;
        }
        inline void set_cost(double c) {
            m_cost = c;
        }
        inline void set_leaf(bool f) {
            m_leaf = f;
        }
        inline bool get_leaf() {
            return m_leaf;
        }
        inline void set_parent(search_node* p) {
            m_parent = p;
        }
        inline size_t num_children() {
            return m_children.size();
        }
        inline search_node* get_parent() {
            return m_parent;
        }
        inline search_node* get_child(size_t v) {
            assert (v >= 0 && v < m_children.size());
            return m_children[v];
        }
        inline void add_children(std::vector<search_node*>& chi) {
            m_children = chi;
        }
        inline void add_child(search_node* c) {
            m_children.push_back(c);
        }
        inline std::vector<search_node*>& get_children() {
            return m_children;
        }
        inline void remove_child(search_node* c) {
            std::vector<search_node*>::iterator it = std::find(m_children.begin(), m_children.end(), c);
            if (it != m_children.end()) {
                m_children.erase(it);
            }
        }
        inline std::vector<double>& get_cache() {
            return m_cache;
        }
        inline void set_cache(std::vector<double>& cache) {
            m_cache = cache;
        }
        inline static bool heur_less(const search_node* a, const search_node* b) {
            return a->get_heur() < b->get_heur();
        }
        inline static bool heur_greater(const search_node* a, const search_node* b) {
            return a->get_heur() > b->get_heur();

        }
        inline std::map<size_t, size_t> get_path_assignment() {
            std::map<size_t, size_t> assgn;
            return assgn;
        }
        inline std::vector<int>& get_assignment() {
            return m_assignment;
        }
        inline void set_assignment(std::vector<int>& assgn) {
            m_assignment = assgn;
        }
        inline void clear_assignment() {
            m_assignment.clear();
        }
        inline bool is_cachable() {
            return m_cachable;
        }
        inline void set_cachable() {
            m_cachable = true;
        }
        inline void set_pruned() {
            m_pruned = true;
        }
        inline bool is_pruned() {
            return m_pruned;
        }
        inline void set_optimal(bool f) {
            m_optimal = f;
        }
        inline bool is_optimal() {
            return m_optimal;
        }
        inline void set_context(std::string context) {
            m_context = context; // the value assignment to the context variables
        }
        inline std::string get_context() {
            return m_context; // the value assignment to the context variables
        }
        inline void set_argmax(int v) {
            m_argmax = v;
        }
        inline int get_argmax() {
            return m_argmax;
        }
        inline std::string to_string() {
            std::ostringstream oss;
            if (m_type == MERLIN_NODE_AND) {
                oss << "AND node: x" << m_variable << "=" << m_value << " [ "
                    << "  w = " << m_weight
                    << ", h = " << m_heur
                    << ", v = " << m_cost
                    << ", c = " << m_children.size()
                    << ", s = " << m_subsolved << " ]"
                    ;
            } else {
                oss << "OR node: x" << m_variable << " [ "
                    << "  h = " << m_heur
                    << ", v = " << m_cost
                    << ", c = " << m_children.size() 
                    << ", a = " << m_argmax << " ]"
                    ;
            
            }

            return oss.str();
        }
    };


} // end namespace

#endif



