/*
 * ipe2u.cpp
 *
 *  Created on: Oct 19, 2020
 *      Author: radu
 *
 * Copyright (c) 2020, International Business Machines Corporation
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


#include "ipe2u.h"

namespace merlin {

// Initialize IPE2U
void ipe2u::init() {

	// Prologue
	std::cout << "[IPE2U] Begin initialization ..." << std::endl;
	std::cout << "[IPE2U] Find minfill elimination order" << std::endl;
	rand_seed(m_seed); // set the random number generator seed
	if (m_order.size() == 0) { // if we need to construct an elimination ordering
		m_order = order();
	}

	// Find minfill elimination order
	std::cout << "[IPE2U] Minfill order: ";
	std::copy(m_order.begin(), m_order.end(), std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;

	// Calculate the induced width of the elimination ordering
	size_t wstar = induced_width(m_order);
	std::cout << "[IPE2U] Induced width: " << wstar << std::endl;

	// Initialize variable processing order (schedule) as a topological order of the graph
	m_schedule = topological_sort();
	std::cout << "[IPE2U] Schedule: ";
	std::copy(m_schedule.begin(), m_schedule.end(), std::ostream_iterator<size_t>(std::cout, " "));
	std::cout << std::endl;

	// Initialize the Pi- and Lambda-messages for each variable
	for (size_t v = 0; v < nvar(); ++v) {
		//variable_set vs(var(v));
		m_pi.push_back(interval(interval::value(1.0, 1.0)));
		m_lambda.push_back(interval(interval::value(1.0, 1.0)));
	}

	// Initialize the parent-to-child and child-to-parent messages
	m_incoming.resize(nvar());
	m_outgoing.resize(nvar());
	std::vector<directed_edge>::const_iterator vi = m_edges.begin();
	for (; vi != m_edges.end(); ++vi) {
		const directed_edge& e = (*vi);
		variable p = var(e.first), c = var(e.second);
		ipe2u::message m(p, c);
		m_messages.push_back(m);
		size_t mi = m_messages.size() - 1; // get the index of the last message
		m_outgoing[e.first].push_back(mi);
		m_incoming[e.second].push_back(mi);
	}

	// Init evidence
	if (m_evidence.empty() == false) {
		size_t d = nvar();
		std::cout << "[IPE2U] Evidence:";
		std::map<size_t, size_t>::iterator ei = m_evidence.begin();
		for (; ei != m_evidence.end(); ++ei) {
			size_t x = ei->first, xi = ei->second;
			std::cout << " x" << x << "=" << xi;
			variable p = var(x);
			variable dummy(d++, 1);
			ipe2u::message m(p, dummy);
			m.lambda = (xi == 0) ? interval::value(0, 0) : interval::value(infty(), infty());
			m.vacuous = true;

			m_messages.push_back(m);
			size_t mi = m_messages.size() - 1; // get the index of the last message
			m_outgoing[x].push_back(mi);
		}
		std::cout << std::endl;
	}

	// Output the incoming and outgoing messages for each variable
	if (m_verbose > 0) {
		std::cout << "[DEBUG] Incoming and outgoing edges:" << std::endl;
		for (size_t v = 0; v < nvar(); ++v) {
			std::cout << " x " << v << ": i(";
			for (size_t i = 0; i < m_incoming[v].size(); ++i) {
				size_t mindex = m_incoming[v][i];
				message& m = m_messages[mindex];
				std::cout << m.parent;
				if (i < m_incoming[v].size() - 1) std::cout << " ";
			}
			std::cout << ")  o(";
			for (size_t j = 0; j < m_outgoing[v].size(); ++j) {
				size_t mindex = m_outgoing[v][j];
				message& m = m_messages[mindex];
				std::cout << m.child;
				if (j < m_outgoing[v].size()) std::cout << " ";
			}
			std::cout << ")" << std::endl;
		}
		std::cout << "[DEBUG] Initial PI messages:" << std::endl;
		for (size_t i = 0; i < m_pi.size(); ++i) {
			std::cout << "x" << i << ": " << m_pi[i] << std::endl;
		}
		std::cout << "[DEBUG] Initial LAMBDA messages:" << std::endl;
		for (size_t i = 0; i < m_lambda.size(); ++i) {
			std::cout << "x" << i << ": " << m_lambda[i] << std::endl;
		}
		std::cout << "[DEBUG] Initial messages:" << std::endl;
		for (size_t i = 0; i < m_messages.size(); ++i) {
			std::cout << i << ": (" << m_messages[i].parent << "->"
				<< m_messages[i].child << "), pi=" << m_messages[i].pi
				<< ", lambda=" << m_messages[i].lambda << std::endl; 
		}
		std::cout << std::endl;
	}

	// Initialize the average change in messages (lower and upper bounds)
	m_delta = interval::value(0.0, 0.0);

	// Calculate total initialization time
	double elapsed = (timeSystem() - m_start_time);
	std::cout << "[IPE2U] Finished initialization in " << elapsed << " seconds" << std::endl;
}

// Update the marginals (after belief propagation)
void ipe2u::update_beliefs() {
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

		if (m_verbose > 0) {
			std::cout << "[DEBUG] Belief for x" << x << ": " << bel << std::endl;
		}

		m_beliefs[x] = bel;
	}
}

// Compute the pi(x) message (x: current) -- update the m_pi[x] factor
void ipe2u::pi(variable x) {

	size_t v = x.label();
	std::vector<variable> parents;
	for (size_t i = 0; i < m_incoming[v].size(); ++i) {
		size_t mi = m_incoming[v][i];
		parents.push_back(m_messages[mi].parent);
	}

	// Check for root node
	if (parents.empty()) { // this is a root node
		m_pi[v][0] = m_factors[v].get(1);
	} else { // this is a node with parents

		double lb = infty(), ub = -infty();
		size_t n = parents.size();

		// Enumerate all possible extreme points of the parent intervals
		std::vector<int> extremes(n, 0);
		extremes[n - 1] = -1;
		int i;
		while (true) {
			for (i = n - 1; i >= 0; --i) {
				size_t last = parents[i].states() - 1;
				if (extremes[i] < (int)last) break;
				extremes[i] = 0;
			}

			if (i < 0) break; // done
			++extremes[i];

			// Complete configuration of the extreme points
			double sum_lb = 0.0, sum_ub = 0.0;
			std::vector<int> values(n, 0);
			values[n - 1] = -1;
			int j;
			while (true) {
				for (j = n - 1; j >= 0; --j) {
					size_t last = parents[j].states() - 1;
					if (values[j] < (int)last) break;
					values[j] = 0;
				}

				if (j < 0) break; // done
				++values[j];

				// We have a complete configuration of the parents and child
				std::map<size_t, size_t> config;
				config[x.label()] = 1; // child's value=1
				for (size_t vi = 0; vi != values.size(); ++vi) {
					config[parents[vi].label()] = values[vi]; // parents values 
				}

				// Get P(X|U1...Un)
				interval::value val = m_factors[v].get_value(config);
				double plb = val.first; // lower bound
				double pub = val.second; // upper bound

				// Get the parent pi-messages
				for (size_t ui = 0; ui < n; ++ui) {
					size_t mi = m_incoming[v][ui];
					ipe2u::message& m = m_messages[mi]; // 
					size_t s = values[ui];
					double pi = 1.0;
					interval::value interval = m.pi[0];
					if (extremes[ui] == 0) { // low
						pi = (s == 1) ? interval.first : (1.0 - interval.first);
					} else { // high
						pi = (s == 1) ? interval.second : (1.0 - interval.second);
					}

					plb *= pi;
					pub *= pi;
				}

				sum_lb += plb;
				sum_ub += pub;
			}

			// minimize/maximize
			lb = std::min(lb, sum_lb);
			ub = std::max(ub, sum_ub);
		}

		assert(lb <= ub); // safety checks
		m_delta.first += fabs(m_pi[v][0].first - lb);
		m_delta.second += fabs(m_pi[v][0].second - ub);
		m_pi[v].set(0, interval::value(lb, ub));
	}

	if (m_verbose > 0) {
		std::cout << "[DEBUG] Computed PI(x) for x" << x << ": " << m_pi[v] << std::endl;
	}
}

// Compute the lambda(x) message (x: current) -- update the m_lambda[x] factor
void ipe2u::lambda(variable x) {

	size_t v = x.label();
	double lb = 1.0, ub = 1.0;
	for (size_t j = 0; j < m_outgoing[v].size(); ++j) {
		size_t mi = m_outgoing[v][j];
		ipe2u::message& m = m_messages[mi];
		lb *= m.lambda[0].first;
		ub *= m.lambda[0].second;
	}

	m_delta.first += (isfinite(lb) && isfinite(m_lambda[v][0].first)) ? fabs(m_lambda[v][0].first - lb) : 0.0;
	m_delta.second += (isfinite(ub) && isfinite(m_lambda[v][0].second)) ? fabs(m_lambda[v][0].second - ub) : 0.0;
	m_lambda[v].set(0, interval::value(lb, ub));

	if (m_verbose > 0) {
		std::cout << "[DEBUG] Computed LAMBDA(x) for x" << x << ": " << m_lambda[v] << std::endl;
	}
}

// Compute the pi(x-child) message (x: current, y: child) -- factor(x)
void ipe2u::pi(variable x, variable y, ipe2u::message& m) {
	
	if (m.evidence || m.vacuous) return; // nothing to do for dummy messages (evidence or vacuous)

	interval& f = m.pi; // pi-message to be updated
	size_t v = x.label();
	double lb = 1.0, ub = 1.0;
	double plb = 1.0, pub = 1.0;
	interval::value Px = m_pi[v].get(0);
	for (size_t j = 0; j < m_outgoing[v].size(); ++j) {
		size_t mi = m_outgoing[v][j];
		ipe2u::message& t = m_messages[mi];
		interval::value L = t.lambda.get(0);
		if (t.child != y) {
			plb *= L.first; // lower bound
			pub *= L.second; // upper bound
		}
	}

	// Compute the lower and upper bounds of the pi-message
	lb = 1.0 / (1.0 + (1.0/Px.first - 1.0) * (1.0/plb));
	ub = 1.0 / (1.0 + (1.0/Px.second - 1.0) * (1.0/pub));

	assert(lb <= ub);
	m_delta.first += fabs(f[0].first - lb);
	m_delta.second += fabs(f[0].second - ub);
	f[0] = interval::value(lb, ub);

	if (m_verbose > 0) {
		std::cout << "[DEBUG] Computed PI(x" << x << ",x" << y << ") for parent x" << x << " to child x" << y << ": " << f << std::endl;
	}
}

double ipe2u::hi(variable x, variable u, size_t ui, bool low, 
	std::vector<variable>& parents, std::vector<int>& extremes) {
	
	size_t v = x.label();
	double result = 0.0;

	// Handle the special case: no parents
	if (parents.empty()) {
		std::map<size_t, size_t> config;
		config[x.label()] = 1; // child's value=1
		config[u.label()] = ui;
		double prod = 1.0;
		if (low == true) {
			prod *= m_factors[v].get_value(config).first; // lower bound
		} else {
			prod *= m_factors[v].get_value(config).second; // upper bound
		}

		result += prod;
		return result;
	}

	size_t n = parents.size();
	std::vector<int> values(n, 0);
	values[n - 1] = -1;
	int j;
	while (true) {
		for (j = n - 1; j >= 0; --j) {
			size_t last = parents[j].states() - 1;
			if (values[j] < (int)last) break;
			values[j] = 0;
		}

		if (j < 0) break; // done
		++values[j];

		// We have a complete configuration of the parents
		std::map<size_t, size_t> config;
		for (size_t vi = 0; vi != values.size(); ++vi) {
			config[parents[vi].label()] = values[vi]; // parents values 
		}

		config[x.label()] = 1; // child's value=1
		config[u.label()] = ui;
		double prod = 1.0;
		if (low == true) {
			prod *= m_factors[v].get_value(config).first; // lower bound
		} else {
			prod *= m_factors[v].get_value(config).second; // upper bound
		}

		for (size_t vi = 0; vi < n; ++vi) {
			size_t mi = m_incoming[v][vi];
			ipe2u::message& m = m_messages[mi]; // 
			size_t s = values[vi];
			double pi = 1.0;
			interval::value interval = m.pi[0];
			if (extremes[vi] == 0) { // low
				pi = (s == 1) ? interval.first : (1.0 - interval.first);
			} else { // high
				pi = (s == 1) ? interval.second : (1.0 - interval.second);
			}
			prod *= pi;
		}

		result += prod;
	} // end while

	return result;
}

double ipe2u::gi1(variable x, variable u, double L, 
	std::vector<variable>& parents, std::vector<int>& extremes) {
	
	double a = hi(x, u, 1, false, parents, extremes);
	double b = hi(x, u, 0, true, parents, extremes);

	if (isinfty(L)) {
		return (b==0.0) ? infty() : (a/b);
	} else {
		return ((L - 1.0)*a + 1.0) / ((L - 1.0)*b + 1.0);
	}
}

double ipe2u::gi2(variable x, variable u, double L, 
	std::vector<variable>& parents, std::vector<int>& extremes) {
	
	double a = hi(x, u, 1, true, parents, extremes);
	double b = hi(x, u, 0, false, parents, extremes);

	if (isinfty(L)) {
		return (b==0.0) ? infty() : (a/b);
	} else {
		return ((L - 1.0)*a + 1.0) / ((L - 1.0)*b + 1.0);
	}
}

// Compute the lambda(x-parent) message (x: current, u: parent) -- scalar
void ipe2u::lambda(variable x, variable u, message& m) {
	
	if (m.vacuous) return; // do nothing for vacuous messages

	interval& f = m.lambda; // lambda-message to be updated
	size_t v = x.label();
	std::vector<variable> parents; // other than 'u'
	for (size_t i = 0; i < m_incoming[v].size(); ++i) {
		size_t mi = m_incoming[v][i];
		variable par = m_messages[mi].parent;
		if (par != u) {
			parents.push_back(par);
		}
	}

	std::vector<double> Lx;
	Lx.push_back(m_lambda[v][0].first);
	Lx.push_back(m_lambda[v][0].second);
	double lb = infty(), ub = -infty(); // min/max
	for (size_t l = 0; l < Lx.size(); ++l) {
		double L = Lx[l];

		if (parents.empty()) {
			std::vector<int> extremes; // an empty vector
			double gi_lb = (L <= 1) ? gi1(x, u, L, parents, extremes) 
				: gi2(x, u, L, parents, extremes);
			double gi_ub = (L <= 1) ? gi2(x, u, L, parents, extremes) 
				: gi1(x, u, L, parents, extremes);

			lb = std::min(lb, gi_lb);
			ub = std::max(ub, gi_ub);
		} else {
			// Enumerate all possible extreme points of the parent intervals
			size_t n = parents.size();
			std::vector<int> extremes(n, 0);
			extremes[n - 1] = -1;
			int i;
			while (true) {
				for (i = n - 1; i >= 0; --i) {
					size_t last = parents[i].states() - 1;
					if (extremes[i] < (int)last) break;
					extremes[i] = 0;
				}

				if (i < 0) break; // done
				++extremes[i];

				// Complete configuration of the extreme points
				double gi_lb = (L <= 1) ? gi1(x, u, L, parents, extremes) 
					: gi2(x, u, L, parents, extremes);
				double gi_ub = (L <= 1) ? gi2(x, u, L, parents, extremes) 
					: gi1(x, u, L, parents, extremes);

				lb = std::min(lb, gi_lb);
				ub = std::max(ub, gi_ub);
			} // end while
		}
	}

	m_delta.first += (isfinite(lb) && isfinite(f[0].first)) ? fabs(f[0].first - lb) : 0.0;
	m_delta.second += (isfinite(ub) && isfinite(f[0].second)) ? fabs(f[0].second - ub) : 0.0;
	f[0] = interval::value(lb, ub);
	if (m_verbose > 0) {
		std::cout << "[DEBUG] Computed LAMBDA(x" << x << ",x" << u << ") from child x" << x << " to parent x" << u << ": " << f << std::endl;
 	}
}

// Message propagation (2U) on a polytree
void ipe2u::propagate() {
	// Message-passing for a number of iterations
	double num_messages = 2.0 * (double)nvar() + 2.0 * edges().size();
	std::cout << "[IPE2U] Messages per iteration: " << num_messages << std::endl;
	std::cout << "[IPE2U] Begin message-passing ..." << std::endl;
	for (size_t iter = 1; iter <= 100; ++iter) {
		m_delta = interval::value(0.0, 0.0);
		for (size_t k = 0; k < m_schedule.size(); ++k) {
			size_t v = m_schedule[k];
			variable x = var(v);

			// compute pi(x) based on pi(parent-x)
			pi(x);

			// compute lambda(x) based on lambda(child-x)
			lambda(x);

			// compute pi(x-child) based on pi(x) and lambda(child-x)
			for (size_t j = 0; j < m_outgoing[v].size(); ++j) {
				ipe2u::message& m = m_messages[m_outgoing[v][j]];
				variable yj = m.child;
				pi(x, yj, m);
			}

			// compute lambda(x-parent) based on lambda(x) and pi(parent-x)
			for (size_t i = 0; i < m_incoming[v].size(); ++i) {
				ipe2u::message& m = m_messages[m_incoming[v][i]];
				variable ui = m.parent;
				lambda(x, ui, m);
			}
		}

		// Keep track of the progress
		m_delta.first /= num_messages;
		m_delta.second /= num_messages;
		std::cout << " iteration " << iter << " ... [" << m_delta.first 
			<< "," << m_delta.second << "]" << std::endl;

		// Check for eary convergence
		if (m_delta.first <= m_threshold && m_delta.second <= m_threshold) {
			std::cout << "[IPE2U] Converged after " << iter << " iterations" << std::endl;
			break;
		}
	}

	std::cout << "[IPE2U] Finished message-passing in " << (timeSystem() - m_start_time)
		<< " seconds" << std::endl;

}

void ipe2u::set_vacuous_messages(my_set<directed_edge>& cutset) {
	for (my_set<directed_edge>::const_iterator i = cutset.begin(); 
		i != cutset.end(); ++i) {
		size_t u = i->first;
		size_t v = i->second;
		assert(m_messages[i->idx].parent == u);
		assert(m_messages[i->idx].child == v);
		m_messages[i->idx].vacuous = true;
		m_messages[i->idx].pi = interval(interval::value(0.0, 1.0));
		m_messages[i->idx].lambda = interval(interval::value(0.0, infty()));
	}
}

void ipe2u::intersect_beliefs() {
	std::vector<interval> beliefs(nvar());
	for (size_t x = 0; x < nvar(); ++x) { // for each variable
		double lb = -infty(), ub = infty();
		for (size_t i = 0; i < m_depot.size(); ++i) {
			interval& bel = m_depot[i][x];
			lb = std::max(lb, bel[1].first);
			ub = std::min(ub, bel[1].second);
		}

		beliefs[x] = interval(variable_set(var(x)));
		beliefs[x][0] = interval::value(1.0-ub, 1.0-lb);
		beliefs[x][1] = interval::value(lb, ub);
	}

	// Update the final beliefs
	m_beliefs = beliefs;
}

/// Run the IPE2U algorithm.
void ipe2u::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();
	std::vector<message> original = m_messages; // save the initial messages

	// Iterating over loop cutsets
	std::cout << "[IPE2U] Begin iterating over loop cutsets ..." << std::endl;
	for (size_t iter = 1; iter <= m_iterations; ++iter) {

		// Find a new loop cutset
		my_set<directed_edge> cutset = find_loop_cutset();
		std::cout << " iteration " << iter << " loop cutset size is " << cutset.size() << std::endl;
		
		// Set the vacuous messages
		m_messages = original;
		set_vacuous_messages(cutset);

		// 2U message propagation
		propagate();

		// Update the beliefs
		update_beliefs();
		std::cout << "[IPE2U] Updated the beliefs: ";
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
		m_depot.push_back(m_beliefs); // store the current marginals (for final intersection)
		m_beliefs.clear();
	}

	intersect_beliefs();
	std::cout << "[IPE2U] Finished in " << (timeSystem() - m_start_time)
		<< " seconds" << std::endl;
	
	// Output the marginals 
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

/// Write the solution to the output stream
void ipe2u::write_solution(std::ostream& out, int output_format) {

	if (output_format == MERLIN_OUTPUT_JSON) {
		out << "{";
		out << " \"algorithm\" : \"loopy2u\", ";
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

