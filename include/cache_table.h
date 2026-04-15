/*
 * cache_table.h
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

/// \file cache_table.h
/// \brief Cache table
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_CACHE_TABLE_H_
#define IBM_LOOPY_CACHE_TABLE_H_

#include "base.h"

typedef std::vector<size_t> context_t;
typedef std::pair<double, std::vector<int>> mapped_value;
typedef std::unordered_map<std::string, mapped_value> context_map;

namespace merlin {

///
/// @brief The cache table structure
///
class cache_table {
public:
    cache_table(size_t var) : m_variable(var) {};
    ~cache_table() {};

    inline void write(std::string ctxt, double v, std::vector<int>& sol) {
        m_table.insert(std::make_pair(ctxt, std::make_pair(v, sol)));
    }
	
    inline mapped_value read(std::string ctxt) const {
        context_map::const_iterator ci = m_table.find(ctxt);
        if (ci == m_table.end()) {
            throw 1; // not found
        } else {
            return ci->second;
        }
    }

protected:
    /// Members

    size_t m_variable;                             ///< The variable
    context_map m_table;                           ///< Cache table
};

} // end namespace

#endif

