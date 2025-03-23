/*
 * bound_propagator.h
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

#ifndef IBM_LOOPY_BOUND_PROPAGATOR_H_
#define IBM_LOOPY_BOUND_PROPAGATOR_H_

#include "base.h"
#include "search_node.h"
#include "search_space.h"
#include "pseudotree.h"
namespace merlin {

///
/// @brief Bound propagator for branch and bound search
///
class bound_propagator {
public:
    bound_propagator();
    ~bound_propagator();

    search_node* propagate(search_node* n, bool solution, search_node* upper_limit = NULL);
    void propagate_tuple(search_node* start, search_node* end);
    void update_solution(double timestamp, double cost, std::vector<int>& sol, std::pair<size_t, size_t> num_nodes);

    inline double get_best_cost() {
        return m_best_cost;
    };
    inline std::vector<int> get_best_config() {
        return m_best_config;
    }
    inline void init(double start_time, pseudotree* pt, search_space* s, bool caching = false) {
        m_start_time = start_time;
        m_pseudotree = pt;
        m_space = s;
        m_caching = caching;
        m_best_cost = NAN; // NaN
    }

protected:
    /// Members

    bool m_caching;                                 ///< Enable caching
    std::vector<std::vector<int>> m_solutions;      ///< Solutions
    double m_best_cost;                             ///< Best solution cost so far
    std::vector<int> m_best_config;                 ///< Best solution so far
    pseudotree* m_pseudotree;                       ///< Pseudo tree
    double m_start_time;                            ///< Start time
    search_space* m_space;                          ///< Search space
};

} // end namespace

#endif

