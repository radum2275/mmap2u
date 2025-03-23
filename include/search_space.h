/*
 * search_space.h
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

/// \file search_space.h
/// \brief Search space
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_SEARCH_SPACE_H_
#define IBM_LOOPY_SEARCH_SPACE_H_

#include "base.h"
#include "search_node.h"
#include "cache_table.h"

namespace merlin {

///
/// @brief The search space structure
///
class search_space {
public:

    /// @brief Constructor
    search_space() {
        m_root = NULL;
        m_num_vars = 0;
        m_num_nodes = std::make_pair(0, 0);
    };

    /// @brief Destructor
    ~search_space() {};

    /// @brief Set the root of the search space.
    /// @param n is a search node
    inline void set_root(search_node* n) {
        assert(n != NULL);
        m_root = n;
    }

    /// @brief Get the root of the search space.
    /// @return a pointer to the root node
    inline search_node* get_root() {
        return m_root;
    }

    /// @brief Get the number of nodes expanded.
    /// @return a pair of two numbers representing the AND and OR nodes.
    inline std::pair<size_t, size_t> get_num_nodes() {
        return m_num_nodes;
    }

    /// @brief Add a node to the search space (count).
    /// @param type is the node type (AND/OR)
    inline void add_node(size_t type) {
        if (type == MERLIN_NODE_AND) {
            m_num_nodes.first += 1;
        } else {
            m_num_nodes.second += 1;
        }
    }

    /// @brief Initialize the search space.
    /// @param num_vars is the number of variables.
    inline void init(size_t num_vars) {
        assert(num_vars > 0);
        m_num_vars = num_vars;
        for (size_t i = 0; i < num_vars; ++i) {
            m_cache.push_back(std::make_unique<cache_table>(i));
        }
    }

    /// @brief Write a node cost assignment pair to a cache table.
    /// @param var is the variable corresponding to the OR node
    /// @param ctxt is the context of the variable
    /// @param v is the cost of the subproblem rooted an the node
    /// @param sol is the optimal assignment of the subproblem rooted at the node 
    inline void write(size_t var, std::string ctxt, double v, std::vector<int>& sol) {
        assert(var >= 0 && var < m_cache.size());
        m_cache[var]->write(ctxt, v, sol);
    }
	
    /// @brief Read a node cost assignment pair from a cache table.
    /// @param var is the variable corresponding to the OR node
    /// @param ctxt is the context of the variable
    /// @return a pair representing the subproblem cost and optimal assignment.
    inline mapped_value read(size_t var, std::string ctxt) const {
        assert(var >= 0 && var < m_cache.size());
        return m_cache[var]->read(ctxt); // throws int exception if entry not found
    }

protected:
    /// Members

    search_node* m_root;                                ///< Root
    size_t m_num_vars;                                  ///< Number of variables
    std::pair<size_t, size_t> m_num_nodes;              ///< Number of nodes (AND, OR)
    std::vector<std::unique_ptr<cache_table>> m_cache;  ///< Cache
};


} // end namespace

#endif

