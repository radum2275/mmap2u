/*
 * util.cpp
 *
 *  Created on: 17 Sep 2018
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

/// \file util.cpp
/// \brief Various utilities
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#include "utils.h"

namespace merlin {

// check if A dominates B (set B is included in set A)
bool dominates(std::set<size_t>& A, std::set<size_t>& B) {
	std::set<size_t> temp;
	std::set_intersection(A.begin(), A.end(),
			B.begin(), B.end(),
			std::inserter(temp, temp.begin()));

	return (temp == B);
}

void rand_seed() {
	srand(time(0));
}
void rand_seed(size_t s) {
	srand(s);
}
// randu returns a random number in [0..1]
double randu() {
	return rand() / double(RAND_MAX);
}
// randi returns a random integer in 0..imax-1
int randi(int imax) {
	assert(imax > 0);
	imax--;
	if (imax == 0)
		return 0;
	int guard = (int) (randu() * imax) + 1;
	return (guard > imax) ? imax : guard;
}
// randi returns a random integer in 0..imax-1
int randi2(int imax) {
	assert(imax > 0);
	int guard = (int) (randu() * imax);
	return (guard >= imax) ? imax-1 : guard;
}
double randn() {  // Marsaglia polar method
	double u, v, s;
	u = 2 * randu() - 1;
	v = 2 * randu() - 1;
	s = u * u + v * v;
	return u * std::sqrt(-2 * std::log(s) / s);
} 

} // end namespace
