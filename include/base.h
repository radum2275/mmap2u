/*
 * base.h
 *
 *  Created on: Feb 8, 2013
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

/// \file base.h
/// \brief Global definitions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

// Software version
#define VERSIONINFO "loopy 1.0.0"
#define COPYRIGHT "(c) Copyright IBM Corp. 2020\nAll Rights Reserved"

#ifndef IBM_MERLIN_BASE_H_
#define IBM_MERLIN_BASE_H_

// Debugging purposes
//#define MERLIN_DEBUG

// Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/timeb.h>

// STL kernel
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <deque>
#include <list>
#include <queue>
#include <set>
#include <stack>
#include <exception>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>
#include <numeric>
#include <cmath>

#define	MERLIN_EXIT_FAILURE	1
#define	MERLIN_EXIT_SUCCESS	0

/// Miscelaneous constants
#define MERLIN_PRECISION 	6			///< Precision used for displaying doubles (default 6)
#define MERLIN_EPSILON 		1e-6		///< Small epsilon value to control determinism
#define MERLIN_UNKNOWN		-1			///< Unknown value


///
/// Probabilistic inference algorithms.
///
#define MERLIN_ALGO_L2U 	    1000		///< Loopy 2U
#define MERLIN_ALGO_IPE2U       1001        ///< Iterated Partial Evaluation 2U
#define MERLIN_ALGO_SV2U        1002        ///< Structural Variational 2U
#define MERLIN_ALGO_CONVERT     1003        ///< Convert to interval credal net
#define MERLIN_ALGO_GENERATOR   1004        ///< Random problem generator
#define MERLIN_ALGO_MMAP_HILL   1005        ///< Stochastic Hill-climbing for MMAP
#define MERLIN_ALGO_MMAP_TABOO  1006        ///< Taboo search for MMAP
#define MERLIN_ALGO_MMAP_SA     1007        ///< Simulated annealing for MMAP
#define MERLIN_ALGO_MMAP_GLS    1008        ///< Guided local search for MMAP

///
/// Probabilistic inference tasks.
///
#define MERLIN_TASK_MAR		20			///< Posterior marginals (given evidence)
#define MERLIN_TASK_CONV    30          ///< Convert/translate 
#define MERLIN_TASK_MMAP    40          ///< Marginal MAP (given evidence)
#define MERLIN_TASK_GEN     50          ///< Problem generator

///
/// Credal MMAP type.
///
#define MERLIN_MMAP_MAXIMIN     0
#define MERLIN_MMAP_MAXIMAX     1
#define MERLIN_MMAP_INTERVAL    2

///
/// Input graphical models.
///
#define MERLIN_INPUT_INTERVAL	1		///< UAI Interval credal nets
#define MERLIN_INPUT_BAYES      2       ///< UAI Bayes net
#define MERLIN_INPUT_MARKOV     3       ///< UAI Markov net

///
/// Output format
///
#define MERLIN_OUTPUT_UAI	10			///< UAI output format (default)
#define MERLIN_OUTPUT_JSON	11			///< JSON output format

#endif /* IBM_MERLIN_BASE_H_ */
