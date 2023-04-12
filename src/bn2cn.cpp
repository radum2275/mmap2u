/*
 * bn2cn.cpp
 *
 *  Created on: 21 Apr 2021
 *      Author: radu
 *
 * Copyright (c) 2021, International Business Machines Corporation
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


#include "bn2cn.h"

namespace merlin {

// Initialize Loopy2U
void bn2cn::init() {

	// Prologue
	std::cout << "[BN2CN] Begin initialization ..." << std::endl;
	std::cout << "[BN2CN] Random generator seed: " << m_seed << std::endl;

	// Calculate total initialization time
	double elapsed = (timeSystem() - m_start_time);
	std::cout << "[BN2CN] Finished initialization in " << elapsed << " seconds" << std::endl;
}

/// Run the BN2CN conversion.
void bn2cn::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();

	// Convert the factor values into intervals
	std::cout << "[BN2CN] Converting the factor values to intervals ..." << std::endl;
	for (size_t i = 0; i < m_factors.size(); ++i) {
		interval& fn = m_factors[i];
		for (size_t j = 0; j < fn.numel(); ++j) {
			interval::value val = fn[j];
			double a = val.first;
			double b = val.second;
			a = std::max(0.0, a - m_epsilon);
			b = std::min(1.0, b + m_epsilon);
			fn[j] = interval::value(a, b);
		}
	}

	std::cout << "[BN2CN] Finished conversion in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
}

/// Write the solution to the output stream
void bn2cn::write_solution(std::ostream& os, int output_format) {

	assert(output_format == MERLIN_OUTPUT_UAI);
	std::cout << "[BN2CN] Writing converted credal net to the output file" << std::endl;
	
	// Write the header
	size_t n = nvar();
	os << "INTERVAL" << std::endl;
	os << n << std::endl;
	for (size_t i = 0; i < n; ++i) {
		os << m_dims[i];
		if (i < n - 1) {
			os << " ";
		}
	}
	os << std::endl;

	// Write the factor scopes
	size_t m = m_factors.size();
	std::vector<std::vector<variable> > scopes;
	assert(n == m); // For Bayes nets and credal nets only
	os << m << std::endl;
	for (size_t i = 0; i < m; ++i) {
		variable_set scope = m_factors[i].vars();
		int child = m_factors[i].get_child();
		std::vector<variable> temp;
		for (variable_set::const_iterator vi = scope.begin();
				vi != scope.end(); ++vi) {
			if (vi->label() != (size_t)child) {
				temp.push_back(*vi);
			}
		}

		temp.push_back(var(child)); // child is last variable
		scopes.push_back(temp);

		os << temp.size();
		for (size_t j = 0; j < temp.size(); ++j) {
			os << " " << temp[j].label();
		}

		os << std::endl;
	}

	os << std::endl;

	// Write the factor tables
	for (size_t i = 0; i < m; ++i) {
		const interval& fn = m_factors[i];
		os << fn.numel() << std::endl;
		std::vector<variable> orig_scope;
		for (variable_set::const_iterator vi = fn.vars().begin();
				vi != fn.vars().end(); ++vi) {
			orig_scope.push_back(*vi);
		}
		convert_index cv(scopes[i], false, orig_scope, true);
		for (size_t j = 0; j < fn.numel(); ++j) {
			size_t k = cv.convert(j);
			interval::value val = fn.get(k);
			os << " " << std::setiosflags(std::ios::fixed)
				<< std::setprecision(8) << val.first << " " << val.second << std::endl;
		}

		os << std::endl << std::endl;
	}
}

} // namespace

