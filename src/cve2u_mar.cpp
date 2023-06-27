/*
 * cve2u_mar.cpp
 *
 *  Created on: 04 Apr 2023
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


#include "cve2u_mar.h"

namespace merlin {

// Initialize CVE2U
void cve2u_mar::init() {

	// Prologue
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Begin initialization ..." << std::endl;
		std::cout << "[CVE2U] Random generator seed: " << m_seed << std::endl;
		std::cout << "[CVE2U] Find minfill elimination order ..." << std::endl;
	}
	
	rand_seed(m_seed); // set the random number generator seed
	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = order();
	}

	// Find minfill elimination order
	std::cout << "[CVE2U] Original minfill order: ";
	std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;

	// Set the query variable at the end of the ordering
	assert(m_query.size() == 1);
	std::vector<size_t>::iterator pos = std::find(m_order.begin(), m_order.end(), m_query[0]);
	m_order.erase(pos);
	m_order.push_back(m_query[0]);

	std::cout << "[CVE2U] Revised minfill order: ";
	std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;

	std::cout << "[CVE2U] Query variable: " << m_query[0] << std::endl;

	// Set the variable positions in the ordering
	m_positions.resize(m_order.size());
	for (size_t i = 0; i < m_order.size(); ++i) {
		size_t var = m_order[i];
		m_positions[var] = i;
	}

	// Calculate the induced width of the elimination ordering
	size_t wstar = induced_width(m_order);
	std::cout << "[CVE2U] Induced width: " << wstar << std::endl;

	// Init evidence
	if (m_evidence.empty() == false) {
		std::cout << "[CVE2U] Evidence:";
		std::map<size_t, size_t>::iterator ei = m_evidence.begin();
		for (; ei != m_evidence.end(); ++ei) {
			size_t x = ei->first, xi = ei->second;
			std::cout << " x" << x << "=" << xi;
		}
		std::cout << std::endl;
	} else {
		std::cout << "[CVE2U] Evidence: none" << std::endl;
	}

	// Calculate total initialization time
	double elapsed = (timeSystem() - m_start_time);
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Finished initialization in " << elapsed << " seconds" << std::endl;
	}
}

// Reset the internal state of the solver
void cve2u_mar::reset() {
	
	// Clear the buckets
	m_buckets.clear();

	// Initialize the buckets
	size_t num_vars = nvar();
	m_buckets.resize(num_vars);
    std::vector<bool> used(num_vars, false);
    for (size_t i = 0; i < m_order.size(); ++i) {
        size_t v = m_order[i];
        m_buckets[i].set_variable(v);
        for (size_t j = 0; j < m_factors.size(); ++j) {
            interval& f = m_factors[j];
            int ch = f.get_child();
            if (used[ch] == true) {
                continue;
            } else {
                // check if the current interval factor contains the bucket var
                if (f.vars().contains(var(v))) {
                    used[ch] = true;
                    m_buckets[i].add_potential(f.to_potential());
                }
            }
        }

		// Add evidence potential (if any)
		if (m_evidence.find(v) != m_evidence.end()) {
			size_t val = m_evidence[v];
			potential pot;
			factor f(variable_set(var(v)), 0.0);
			f[val] = 1.0;
			pot.add_p(f);
			pot.add_q(f);
			m_buckets[i].add_potential(pot);
		}

		// Add the query variable potential (z=1 by default)
		if (v == m_query[0]) {
			size_t val = 1;
			potential pot;
			factor f1(variable_set(var(v)), 0.0);
			factor f2(variable_set(var(v)), 0.0);
			f1[0] = 1.0;
			f2[1] = 1.0;
			pot.add_p(f1);
			pot.add_q(f2);
			m_buckets[i].add_potential(pot);
		}
    }
}

// Evaluate a partial variable assignment
std::pair<double, double> cve2u_mar::eval() {
    
	// Reset the solver
	reset();

	// Compute the lower probability of evidence
	double lower = ve(false);

	// Reset the solver
	reset();

	// Compute the upper probability of evidence
	double upper = ve(true);

	// Return the result
	return std::make_pair(lower, upper);
}

// Variable elimination to compute P(e)
double cve2u_mar::ve(bool upper) {

	if (m_verbose > 0) {
		if (upper) {
			std::cout << "[CVE2U] Computing Upper P(x|evidence) ..." << std::endl;
		} else {
			std::cout << "[CVE2U] Computing Lower P(x|evidence) ..." << std::endl;
		}
	}

	size_t num_vars = nvar();

    if (m_verbose > 0) {
        std::cout << "[DEBUG] Bucket structure:" << std::endl;
        for (size_t i = 0; i < m_buckets.size(); ++ i) {
            std::cout << "Bucket [" << m_buckets[i].get_variable() << "]" << std::endl;
            std::vector<potential>& pots = m_buckets[i].potentials(); 
            for (size_t j = 0; j < pots.size(); ++j) {
                std::cout << pots[j] << std::endl;
            } 
        }
    }

	// Eliminate the variables along the minfill ordering
    for (size_t i = 0; i < num_vars; ++i) {
        size_t v = m_order[i];
        variable vx = var(v);

		// Check if the last variable (i.e., query variable)
		if (v == m_query[0]) {
			std::cout << "[CVE] Reached the query (last) variable: " << v << std::endl;
			break; // don't eliminate the query variable
		}

       	std::cout << "[CVE] Eliminating variable: " << v << std::endl;
		if (m_verbose > 0) {
			std::cout << "[DEBUG] Bucket of variable: " << v << std::endl;
            std::vector<potential>& pots = m_buckets[i].potentials(); 
            for (size_t j = 0; j < pots.size(); ++j) {
                std::cout << pots[j] << std::endl;
            } 
		}

        // Combine the potentials in the bucket
        potential comb(1.0);
        std::vector<potential>& pots = m_buckets[i].potentials();
        for (size_t j = 0; j < pots.size(); ++j) {
            comb.multiply2(pots[j]);
        }

        // Eliminate the bucket variable
        potential result = comb;
        result.sum2(vx);

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Result before pruning:" << std::endl;
            std::cout << result << std::endl;
        }

        // Remove dominated vertices
		if (upper == true) {
			result.maximize2();
		} else {
        	result.minimize2();
		}

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Result after pruning:" << std::endl;
            std::cout << result << std::endl;
        }

        // Place new potential in the appropriate bucket
        if (result.isscalar()) {
            m_buckets[num_vars - 1].add_potential(result);
        } else {
            // Find the closest bucket that contains a variable in the potential's scope
            for (size_t j = i + 1; j < num_vars; ++j) {
                int y = m_buckets[j].get_variable();
                variable vy = var(y);
                if (result.vars().contains(vy)) {
                    m_buckets[j].add_potential(result);
                    break;
                }
            }
        }
    }

	// Output the content of the query (last) variable in the ordering
	if (m_verbose > 0) {
		std::cout << "[DEBUG] Bucket of variable: " << m_order[num_vars - 1] << std::endl;
		std::vector<potential>& pots = m_buckets[num_vars - 1].potentials(); 
		for (size_t j = 0; j < pots.size(); ++j) {
			std::cout << pots[j] << std::endl;
		} 
	}

	// After elimination, combine all scalars and potentials of the query variable
    potential r(1.0);
	std::vector<potential>& pots = m_buckets[num_vars - 1].potentials();
	for (size_t j = 0; j < pots.size(); ++j) {
		r.multiply2(pots[j]);
	}

	if (m_verbose > 0) {
		std::cout << "[DEBUG] Final potential:" << std::endl;
		std::cout << r << std::endl;
	}

    // Compute the lower/upper marginal probability of the query variable
	double best_score = (upper ? -infty() : infty());
	variable vx = var(m_order[num_vars - 1]);
	size_t num_elems = r.p().size();
	if (upper == true) {
		for (size_t j = 0; j < num_elems; ++j) {
			factor p = r.get_p(j).sum(variable_set(vx));
			factor q = r.get_q(j).sum(variable_set(vx));
			factor t = q / (p + q);
			if (t[0] > best_score) {
				best_score = t[0];
			}
		}
	} else {
		for (size_t j = 0; j < num_elems; ++j) {
			factor p = r.get_p(j).sum(variable_set(vx));
			factor q = r.get_q(j).sum(variable_set(vx));
			factor t = q / (p + q);
			if (t[0] < best_score) {
				best_score = t[0];
			}
		}
	}

    // Get the best score
	return best_score;
}

/// Run the CVE2U algorithm to compute probability of evidence.
void cve2u_mar::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();

	// Get the lower and upper bounds on P(e)
	reset();
	double lower = ve(false);
	reset();
	double upper = ve(true);

	std::cout << "[CVE2U] Lower P(x|e): " << lower << " (" << std::log10(lower) << ")" << std::endl;
	std::cout << "[CVE2U] Upper P(x|e): " << upper << " (" << std::log10(upper) << ")" << std::endl;
	std::cout << "[CVE2U] Finished in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
}

/// Write the solution to the output stream
void cve2u_mar::write_solution(std::ostream& out, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"cve2u\", ";
		out << " \"task\" : \"MAR\", ";
		out << " \"marginals\" : [ ";

		for (vindex v = 0; v < nvar(); ++v) {
			variable x = var(v);
			out << "{";
			out << " \"variable\" : " << x.label() << ", ";
			out << " \"states\" : " << x.states() << ", ";
			out << " \"probabilities\" : [";
			
			const interval& bel = belief(x);
			for (size_t k = 0; k < x.states(); ++k) {
				out << "{";
				out << " \"state\" : " << k << ", ";
				out << " \"bounds\" : ";
				out << std::fixed << std::setprecision(MERLIN_PRECISION)
					<< "[" << bel[k].first << "," << bel[k].second << "]";
				out << "}";
				if (k < x.states() -1) {
					out << ",";
				}
			}
			
			out << "]} ";
			if (v < nvar() - 1) {
				out << ",";
			}
		}
		out << "]";
		out << "}\n";
	} else if (output_format == MERLIN_OUTPUT_UAI) {
		out << "MAR" << std::endl;
		out << nvar();
		for (vindex v = 0; v < nvar(); ++v) {
			variable x = var(v);
			const interval& bel = belief(x);
			out << " " << x.states();
			for (size_t k = 0; k < x.states(); ++k) {
				out << " " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< "[" << bel[k].first << "," << bel[k].second << "]";
			}
		}
		out << std::endl;
	} else {
		std::string err_msg("[ERROR] Unknown output format.");
		throw std::runtime_error(err_msg);
	}
}

} // namespace

