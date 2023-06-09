/*
 * util.h
 *
 *  Created on: 24 Mar 2015
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

/// \file util.h
/// \brief Various utilities
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_MERLIN_UTIL_H_
#define IBM_MERLIN_UTIL_H_

#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <limits>

#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>

namespace merlin {

inline bool isfinite(double v) {
	return (v <= std::numeric_limits<double>::max() && v >= -std::numeric_limits<double>::max());
}

inline bool isinfty(double v) {
	return !isfinite(v);
}

//inline bool isnan(double v)    { return (v!=v); }

inline double infty() {
	return std::numeric_limits<double>::infinity();
}


// Returns system (wall clock) time in seconds
inline double timeSystem() {
#ifdef WINDOWS
    SYSTEMTIME tbuf;
    GetSystemTime(&tbuf);
    return( (double)(tbuf.wSecond + (double)tbuf.wMilliseconds / 1000.0) );
#else
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return( (double)(tv.tv_sec + (double)tv.tv_usec / 1000000.0) );
#endif
}

inline std::string timestamp() {
	struct tm * dt;
	char buffer[30];
	time_t t = time(NULL);
	dt = localtime(&t);
	strftime(buffer, sizeof(buffer), "[%F %X]", dt);
	return std::string(buffer);
}

inline double timeProcess() {
#ifdef WINDOWS
    throw std::runtime_error("No process time implemented for Windows");
#else
    clock_t tv( clock() );
    return( (double)tv );
#endif
}

//
// String splitting functions
//
inline std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

//
// Random number classes
//
void rand_seed();
void rand_seed(size_t s);

// Returns a random number in [0..1]
double randu();

// Returns a random integer in 0..imax-1
int randi(int imax);

// Returns a random integer in 0..imax-1
int randi2(int imax);

// Marsaglia polar method
double randn();

// check if A dominates B (set B is included in set A)
bool dominates(std::set<size_t>& A, std::set<size_t>& B);

} // namespace



#endif // re-include

