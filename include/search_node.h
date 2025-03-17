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

#define NODE_AND 1
#define NODE_0R 2

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

        std::unique_ptr<search_node> m_parent;
        std::vector<std::unique_ptr<search_node>> m_children;

    public:
        search_node(size_t var, int val, size_t type) 
            : m_variable(var), m_value(val), m_type(type), m_heur(0), m_cost(0), m_parent(nullptr) {};
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
        inline void set_weight(double w) {
            m_weight = w;
        }
        inline void set_depth(size_t d) {
            m_depth = d;
        }
        inline void set_cost(double c) {
            m_cost = c;
        }
        inline void set_parent(std::unique_ptr<search_node> p) {
            m_parent = std::move(p);
        }
        inline size_t num_children() {
            return m_children.size();
        }
        inline search_node& get_parent() {
            return *m_parent;
        }
        inline search_node& get_child(size_t v) {
            assert (v >= 0 && v < m_children.size());
            return *(m_children[v]);
        }
        inline static bool heur_less(const search_node* a, const search_node* b) {
            return a->get_heur() < b->get_heur();
        }
        inline static bool heur_greater(const search_node* a, const search_node* b) {
            return a->get_heur() > b->get_heur();

        }
        void get_path_assignment(std::vector<int>& assignment) {

        }
    };


} // end namespace

#endif



