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

void cve2u::test() {
	size_t n = num_nodes();
	graph g(n); // create the undirected graph
	my_set<directed_edge>::const_iterator ei = m_edges.begin();
	for (; ei != m_edges.end(); ei++) {
		g.add_edge(ei->first, ei->second);
	}

	// Check if graph has cycles
	if (g.is_cyclic()) {
		std::cout << "[DEBUG] Undirected graph is cyclic" << std::endl;
		for (int i = 1; i <= 4; ++i) {
			my_set<directed_edge> cutset = find_loop_cutset();
			std::cout << "[DEBUG] Loop cutset size = " << cutset.size() << std::endl;
			for (my_set<directed_edge>::const_iterator ei = cutset.begin(); ei != cutset.end(); ++ei) {
				std::cout << "[DEBUG] edge (" << ei->first << "," << ei->second << ")" << std::endl;
			}
		}
	} else {
		std::cout << "[DEBUG] Undirected graph is not cyclic (polytree)" << std::endl;
	}
}

// Initialize CVE2U
void cve2u::init() {

	// Prologue
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Begin initialization ..." << std::endl;
		std::cout << "[CVE2U] Random generator seed: " << m_seed << std::endl;
		std::cout << "[CVE2U] Find minfill elimination order" << std::endl;
	}
	rand_seed(m_seed); // set the random number generator seed
	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = order();
	}

	// Find minfill elimination order
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Minfill order: ";
		std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << std::endl;
	}

	// Calculate the induced width of the elimination ordering
	size_t wstar = induced_width(m_order);
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Induced width: " << wstar << std::endl;
	}

	// Init evidence
	if (m_evidence.empty() == false) {
		size_t d = nvar();
		if (m_verbose > 0) {
			std::cout << "[CVE2U] Evidence:";
		}

		std::map<size_t, size_t>::iterator ei = m_evidence.begin();
		for (; ei != m_evidence.end(); ++ei) {
			size_t x = ei->first, xi = ei->second;
			std::cout << " x" << x << "=" << xi;
			variable p = var(x);
			variable dummy(d++, 1);
			message m(p, dummy);
			m.lambda = (xi == 0) ? interval::value(0, 0) : interval::value(infty(), infty());
			m.evidence = true;

			m_messages.push_back(m);
			size_t mi = m_messages.size() - 1; // get the index of the last message
			m_outgoing[x].push_back(mi);
		}
		std::cout << std::endl;
	} else {
		if (m_verbose > 0) {
			std::cout << "[CVE2U] Evidence: none" << std::endl;
		}
	}

	// Calculate total initialization time
	double elapsed = (timeSystem() - m_start_time);
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Finished initialization in " << elapsed << " seconds" << std::endl;
	}
}

// Reset the internal state of the solver
void cve2u::reset() {
	if (m_verbose > 0) {
		std::cout << "[CVE] Reset solver ..." << std::endl;
	}

	m_incoming.clear();
	m_outgoing.clear();
	m_messages.clear();
	m_schedule.clear();
	m_pi.clear();
	m_lambda.clear();
	m_beliefs.clear();
}

// Update the marginals (after belief propagation)
void cve2u::update_beliefs() {
	assert(m_beliefs.empty());
	m_beliefs.resize(nvar());
	for (size_t x = 0; x < nvar(); ++x) {
		variable_set vs(var(x));
		interval bel(vs);
		try { // evidence variable
			size_t k = m_evidence.at(x);
			bel[k] = interval::value(1.0, 1.0);
			bel[1 - k] = interval::value(0.0, 0.0);
		} catch(std::out_of_range e) { // non-evidence variable
			interval::value Lx = m_lambda[x].get(0);
			interval::value Px = m_pi[x].get(0);
			double lb = 1.0 / (1.0 + (1.0 / Px.first - 1.0) * (1.0 / Lx.first));
			double ub = 1.0 / (1.0 + (1.0 / Px.second - 1.0) * (1.0 / Lx.second));
			assert(lb <= ub); // safety checks
			bel[1] = interval::value(lb, ub); // P(x=k|e)
			bel[0] = interval::value(1.0-ub, 1.0-lb); // P(x=0|e)
		}

		if (m_verbose > 1) {
			std::cout << "[DEBUG] Belief for x" << x << ": " << bel << std::endl;
		}

		m_beliefs[x] = bel;
	}

	if (m_verbose > 0) {
		std::cout << "[L2U] Updated the beliefs" << std::endl; 
	}
}

// Compute the lambda(x) message (x: current) -- update the m_lambda[x] factor
void cve2u::lambda(variable x) {

	size_t v = x.label();
	double lb = 1.0, ub = 1.0;
	for (size_t j = 0; j < m_outgoing[v].size(); ++j) {
		size_t mi = m_outgoing[v][j];
		lb *= m_messages[mi].lambda[0].first;
		ub *= m_messages[mi].lambda[0].second;
	}

	m_delta.first += (isfinite(lb) && isfinite(m_lambda[v][0].first)) ? fabs(m_lambda[v][0].first - lb) : 0.0;
	m_delta.second += (isfinite(ub) && isfinite(m_lambda[v][0].second)) ? fabs(m_lambda[v][0].second - ub) : 0.0;
	m_lambda[v].set(0, interval::value(lb, ub));

}


/// Run the Loopy2U algorithm.
void cve2u::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();

	// Eliminate the variables along the minfill ordering

	// Update the beliefs
	update_beliefs();
	if (m_verbose > 0) {
		std::cout << "[CVE2U] Finished in " << (timeSystem() - m_start_time)
			<< " seconds" << std::endl;
	}

	// Output the marginals
	if (m_verbose > 0) {
		std::cout << "MAR" << std::endl;
		std::cout << nvar();
		for (vindex v = 0; v < nvar(); ++v) {
			variable x = var(v);
			const interval& bel = belief(x);
			std::cout << " " << x.states();
			for (size_t k = 0; k < x.states(); ++k) {
				std::cout << " " << std::fixed
					<< std::setprecision(MERLIN_PRECISION)
					<< "[" << bel[k].first << "," << bel[k].second << "]";
			}
		}
	
		std::cout << std::endl;
	}
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

