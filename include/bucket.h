/*
 * bucket.h
 *
 *  Created on: Apr 11, 2023
 *      Author: radu
 *
 * Copyright (c) 2023, International Business Machines Corporation
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


/// \file bucket.h
/// \brief A bucket structure for credal variable elimination
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MMAP2U_BUCKET_H_
#define IBM_MMAP2U_BUCKET_H_

#include "potential.h"

namespace merlin {


///
/// @brief A bucket structure used by Credal Variable Elimination for P(e)
///
class bucket {
public:
    bucket() : m_variable(-1) {};
    ~bucket() {};
    void set_variable(int v) {
        m_variable = v;
    }
    int get_variable() {
        return m_variable;
    }
    void add_potential(const potential& p) {
        m_potentials.push_back(p);
    }
    std::vector<potential>& potentials() {
        return m_potentials;
    }

    std::vector<potential> create_partition(int ibound);

protected:
    int m_variable;
    std::vector<potential> m_potentials;
};

} // end namespace

#endif