/*
 * cve2u.cpp
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


#include "cve2u.h"

namespace merlin {

// Initialize CVE2U
void cve2u::init() {

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
	std::cout << "[CVE2U] Minfill order: ";
	std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;

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
void cve2u::reset() {
	
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
			m_buckets[i].add_potential(pot);
		}
    }
}

// Reset the internal state of the solver
void cve2u::reset(const std::map<size_t, size_t>& config) {
	
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
			m_buckets[i].add_potential(pot);
		}
    }

	// Add the variable assignment as evidence potentials in buckets
	for (std::map<size_t, size_t>::const_iterator mi = config.begin(); mi != config.end(); ++mi) {
		size_t v = mi->first;
		size_t val = mi->second;
		potential pot;
		factor f(variable_set(var(v)), 0.0);
		f[val] = 1.0;
		pot.add_p(f);
		size_t i = m_positions[v];
		m_buckets[i].add_potential(pot);
	}
}

// Evaluate a partial variable assignment
std::pair<double, double> cve2u::eval(const std::map<size_t, size_t> &config) {
    
	// Reset the solver
	reset(config);

	// Compute the lower probability of evidence
	double lower = ve(false);

	// Reset the solver
	reset(config);

	// Compute the upper probability of evidence
	double upper = ve(true);

	// Return the result
	return std::make_pair(lower, upper);
}

// Variable elimination to compute P(e)
double cve2u::ve(bool upper) {

	if (m_verbose > 0) {
		if (upper) {
			std::cout << "[CVE2U] Computing Upper Probability of Evidence ..." << std::endl;
		} else {
			std::cout << "[CVE2U] Computing Lower Probability of Evidence ..." << std::endl;
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
    std::vector<potential> scalars;
    for (size_t i = 0; i < num_vars; ++i) {
        size_t v = m_order[i];
        variable vx = var(v);

		if (m_verbose > 0) {
        	std::cout << "[CVE] Eliminating variable: " << v << std::endl;
		}

        // Combine the potentials in the bucket
        potential comb(1.0);
        std::vector<potential>& pots = m_buckets[i].potentials();
        for (size_t j = 0; j < pots.size(); ++j) {
            comb.multiply(pots[j]);
        }

        // Eliminate the bucket variable
        potential result = comb;
        result.sum(vx);

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Result before pruning:" << std::endl;
            std::cout << result << std::endl;
        }

        // Remove dominated vertices
		if (upper == true) {
			result.maximize();
		} else {
        	result.minimize();
		}

        if (m_verbose > 0) {
            std::cout << "[DEBUG] Result after pruning:" << std::endl;
            std::cout << result << std::endl;
        }

        // Place new potential in the appropriate bucket
        if (result.isscalar()) {
            scalars.push_back(result);
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

    // After elimination, combine all scalars
    potential r(1.0);
    for (size_t i = 0; i < scalars.size(); ++i) {
        r.multiply(scalars[i]);
    }
    
    // Prune dominated scalars
	if (upper == true) {
		r.maximize();
	} else {
    	r.minimize();
	}

    // Check for singleton
    if (r.p().size() > 1) {
        std::cout << "[CVE2U] WARNING: more than one final scalars detected: " << r.p().size() << std::endl; 
    }

    // Get the best score
    double best_score = r.p()[0][0];
	return best_score;
}

/// Run the CVE2U algorithm to compute probability of evidence.
void cve2u::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();

	// Get the lower and upper bounds on P(e)
	reset();
	double lower = ve(false);
	reset();
	double upper = ve(true);

	std::cout << "[CVE2U] Lower P(e): " << lower << " (" << std::log10(lower) << ")" << std::endl;
	std::cout << "[CVE2U] Upper P(e): " << upper << " (" << std::log10(upper) << ")" << std::endl;
	std::cout << "[CVE2U] Finished in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
}

/// Write the solution to the output stream
void cve2u::write_solution(std::ostream& out, int output_format) {

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

