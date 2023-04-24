/*
 * generator.cpp
 *
 *  Created on: 31 Dec 2021
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


#include "generator.h"

namespace merlin {

// Initialize Loopy2U
void generator::init() {

	// Prologue
	std::cout << "[GEN] Begin initialization ..." << std::endl;
	std::cout << "[GEN] Random generator seed: " << m_seed << std::endl;
	std::cout << "[GEN] Graph type: " << m_graph_type << std::endl;
	std::cout << "[GEN] Number of nodes: " << m_num_nodes << std::endl;
	std::cout << "[GEN] Number of parents: " << m_num_parents << std::endl;
	std::cout << "[GEN] KSize: " << m_ksize << std::endl;
	std::cout << "[GEN] Percent: " << m_kpercent << std::endl;
	rand_seed(m_seed);
	m_query_perc = 0;

	// Calculate total initialization time
	double elapsed = (timeSystem() - m_start_time);
	std::cout << "[GEN] Finished initialization in " << elapsed << " seconds" << std::endl;
}

/// Run the generator conversion.
void generator::run() {

	// Start the timer and store it
	m_start_time = timeSystem();

	// Initialize the algorithm
	init();

	// Convert the factor values into intervals
	std::cout << "[GEN] Generating " << m_num_instances << " instances ..." << std::endl;
	for (size_t i = 0; i < m_num_instances; ) {
		cleanup();
		std::ostringstream oss;
		oss << m_graph_type << "_n" << m_num_nodes << "_p" << m_num_parents << "_i" << i << ".uai";
		std::string file_name = oss.str();
		std::cout << "[GEN] Instance: " << file_name << std::endl;
		if (m_graph_type.compare("random") == 0) {
			if (make_random()) {
				std::ofstream os(file_name.c_str());
				write_solution(os, MERLIN_OUTPUT_UAI);
				os.close();
				++i;
			} else {
				std::cout << "[GEN] disconnected!" << std::endl;
				continue;
			}
		} else if (m_graph_type.compare("grid") == 0) {
			make_grid();
			std::ofstream os(file_name.c_str());
			write_solution(os, MERLIN_OUTPUT_UAI);
			os.close();
			++i;


		} else if (m_graph_type.compare("ktree") == 0) {
			make_ktree();
			std::ofstream os(file_name.c_str());
			write_solution(os, MERLIN_OUTPUT_UAI);
			os.close();
			++i;
		} else if (m_graph_type.compare("random_lcn") == 0) {
			if (make_random_lcn3()) {
				std::ofstream os(file_name.c_str());
				write_solution(os, MERLIN_OUTPUT_UAI);
				os.close();
				++i;
			} else {
				std::cout << "[GEN] disconnected!" << std::endl;
				continue;
			}
		} else if (m_graph_type.compare("instance") == 0) {
			// Generating random queries for a given instance
			file_name = m_input_filename;
			std::cout << "[GEN] random queries for instance: " << file_name << std::endl;
			std::ifstream is(file_name);
			read(is);

			size_t num_vars = nvar();
			m_query_perc = m_num_query;
			size_t num_query = (size_t) (num_vars) * (m_num_query / 100.0); // percentage of query vars
			m_num_query = num_query;
			m_num_nodes = num_vars;
			++i;
		}

		// Generate the sample MAP queries
		for (size_t s = 0; s < m_num_samples; ++s) {
			// Select MAP variables at random
			std::set<size_t> query;
			for (size_t q = 0; q < m_num_query; ) {
				size_t v = randi2(m_num_nodes);
				if (query.find(v) == query.end()) {
					query.insert(v);
					++q;
				}
			}

			// Save the query MAP variables in a file
			if (m_graph_type.compare("instance") == 0) {
				std::ostringstream oss2;
				std::string file_name = m_input_filename;
				file_name.erase(file_name.size() - 4); // drop the .uai extension
				oss2 << file_name << ".Q" << m_query_perc << ".s" << s << ".query";
				std::string query_file_name = oss2.str();
				std::ofstream os2(query_file_name.c_str());
				os2 << query.size();
				for (std::set<size_t>::iterator si = query.begin(); si != query.end(); ++si) {
					os2 << " " << *si;
				}
				os2 << std::endl;
				os2.close();

			} else {
				std::ostringstream oss2;
				oss2 << m_graph_type << "_n" << m_num_nodes << "_p" << m_num_parents << "_i" << i-1 << ".Q" << query.size() << ".s" << s << ".query";
				std::string query_file_name = oss2.str();
				std::ofstream os2(query_file_name.c_str());
				os2 << query.size();
				for (std::set<size_t>::iterator si = query.begin(); si != query.end(); ++si) {
					os2 << " " << *si;
				}
				os2 << std::endl;
				os2.close();

				if (m_graph_type.compare("random_lcn") == 0) {
					std::set<size_t> evid;
					for (size_t e = 0; e < m_num_evid;) {
						size_t v = randi2(m_num_nodes);
						if (query.find(v) == query.end() & evid.find(v) == evid.end()) {
							evid.insert(v);
							++e;
						}
					}

					std::ostringstream oss3;
					oss3 << m_graph_type << "_n" << m_num_nodes << "_p" << m_num_parents << "_i" << i-1 << ".E" << evid.size() << ".s" << s << ".evid";
					std::string evid_file_name = oss3.str();
					std::ofstream os3(evid_file_name.c_str());
					os3 << evid.size();
					for (std::set<size_t>::iterator si = evid.begin(); si != evid.end(); ++si) {
						os3 << " " << *si << " 1";
					}
					os3 << std::endl;
					os3.close();
				}
			}
		}
	}

	std::cout << "[GEN] Finished in " << (timeSystem() - m_start_time) << " seconds" << std::endl;
}

void generator::cleanup() {
	m_factors.clear();
	m_vadj.clear();
	m_dims.clear();
}

// Genereate a random network: return true is connected, and false otherwise
bool generator::make_random() {

	// Declare variables.
	size_t n = m_num_nodes;
	size_t p = m_num_parents; // parents per CPT
	size_t r = p + 1; // root CPTs
	size_t c = n - r; // non-root CPTs

	// Create a random ordering of the variables.
	std::vector<size_t> ordering(n);
	std::vector<size_t> position(n);
	for (size_t i = 0; i < n; ++i) {
		ordering[i] = i;
		position[i] = i;
	}

	// Randomly, switch pairs of variables.
	for (size_t i = 0 ; i < n ; ++i) {
		size_t j = randi2(n);
	
		// Switch variable i and j.
		int k = ordering[j];
		ordering[j] = ordering[i];
		ordering[i] = k;

		position[ordering[j]] = j;
		position[ordering[i]] = i;
	}

	// Create the graph structure
	std::vector<std::vector<size_t> > scopes(n);
	std::vector<bool> cpts(n, false);
	size_t count = 0;
	while (count < c) {
		// Pick the child to create P(child|parents).
		int child = ordering[randi2(n - p)];

		// Check if variable was visited
		if (cpts[child]) continue;
		cpts[child] = true;

		// Pick P parents for the child.
		int num_higher_vars = n - position[child] - 1;
		// Notice : number of parents cannot be larger than 'num_higher_vars'.
		size_t num_parents = p;
		if (num_parents > num_higher_vars) 
			num_parents = num_higher_vars;
		
		std::set<size_t> parents;
		for (size_t i = 0; i < num_parents;) {
			size_t parent = ordering[position[child] + 1 + randi2(num_higher_vars)];
			// Check that this parent is not the same as the child.
			if (child == parent) continue;
			// Check other parents.
			if (parents.find(parent) == parents.end()) {
				parents.insert(parent);
				++i;
			}
		}

		// Create a function with child and parents.
		std::vector<size_t> scope(1+num_parents);

		// Set the function's parents.
		int pos = 0;
		std::set<size_t>::iterator it = parents.begin();
		for (; it != parents.end(); ++it) {
			scope[pos++] = (*it);
		}

		// Set the function's child (last position).
		scope[pos] = child;

		// Add the function to the problem instance.
		scopes[child] = scope;
		++count;
	}

	// Create priors (root CPTs).
	for (size_t i = 0; i < n; ++i) {
		if (cpts[i]) continue;
		
		std::vector<size_t> scope(1);
		scope[0] = i;

		scopes[i] = scope;			
	}

	// Create an undirected graph and check if disconnected.
	graph g(n);
	g.init(scopes);
	if (g.is_disconnected()) {
		return false;
	}

	// Create the interval factors 
	for (size_t i = 0; i < scopes.size(); ++i) {
		std::vector<size_t>& scope = scopes[i]; // 'i' is the child var
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		assert(i == (int)scope.back()); // safety check
		f.set_child(scope.back());
		f.fill_random_bayes();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			a = std::max(0.0, a - m_epsilon);
			b = std::min(1.0, b + m_epsilon);
			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

bool generator::make_random_lcn() {

	// Declare variables.
	size_t n = m_num_nodes;
	size_t p = m_num_parents; // max parents per CPT
	size_t r = p + 1; // root CPTs
	size_t c = n - r; // non-root CPTs

	// Create a random ordering of the variables.
	std::vector<size_t> ordering(n);
	std::vector<size_t> position(n);
	for (size_t i = 0; i < n; ++i) {
		ordering[i] = i;
		position[i] = i;
	}

	// Randomly, switch pairs of variables.
	for (size_t i = 0 ; i < n ; ++i) {
		size_t j = randi2(n);
	
		// Switch variable i and j.
		int k = ordering[j];
		ordering[j] = ordering[i];
		ordering[i] = k;

		position[ordering[j]] = j;
		position[ordering[i]] = i;
	}

	// Create the graph structure
	std::vector<std::vector<size_t> > scopes(n);
	std::vector<bool> cpts(n, false);
	size_t count = 0;
	while (count < c) {
		// Pick the child to create P(child|parents).
		int child = ordering[randi2(n - p)];

		// Check if variable was visited
		if (cpts[child]) continue;
		cpts[child] = true;

		// Pick at most P parents for the child (i.e., 1, 2, ... P)
		int num_higher_vars = n - position[child] - 1;
		// Notice : number of parents cannot be larger than 'num_higher_vars'.
		size_t num_parents = randi2(p) + 1;
		if (num_parents > num_higher_vars) 
			num_parents = num_higher_vars;
		
		std::set<size_t> parents;
		for (size_t i = 0; i < num_parents;) {
			size_t parent = ordering[position[child] + 1 + randi2(num_higher_vars)];
			// Check that this parent is not the same as the child.
			if (child == parent) continue;
			// Check other parents.
			if (parents.find(parent) == parents.end()) {
				parents.insert(parent);
				++i;
			}
		}

		// Create a function with child and parents.
		std::vector<size_t> scope(1+num_parents);

		// Set the function's parents.
		int pos = 0;
		std::set<size_t>::iterator it = parents.begin();
		for (; it != parents.end(); ++it) {
			scope[pos++] = (*it);
		}

		// Set the function's child (last position).
		scope[pos] = child;

		// Add the function to the problem instance.
		scopes[child] = scope;
		++count;
	}

	// Create priors (root CPTs).
	for (size_t i = 0; i < n; ++i) {
		if (cpts[i]) continue;
		
		std::vector<size_t> scope(1);
		scope[0] = i;

		scopes[i] = scope;			
	}

	// Create the extra knowledge in the form P(x) or P(x|y) where x,y are two
	// propositional variables. These additional relationships may introduce cycles.
	size_t num_extras = m_num_extras;
	for (size_t e = 0; e < num_extras; ++e) {
		double prob = randu();
		if (prob < 0.0) {
			int x = randi2(n);
			std::vector<size_t> scope(1);
			scope[0] = x;
			scopes.push_back(scope);
		} else {
			int x = randi2(n), y;
			do {
				y = randi2(n);
			} while (y == x);
			std::vector<size_t> scope(2);
			scope[0] = y;
			scope[1] = x;
			scopes.push_back(scope);
		}
	}

	// Create a dummy CPT P(dummy | 1,2....k)
	size_t kk = n-1;
	std::vector<size_t> dummy;
	for (size_t i = 0; i < n-1; ++i) {
		dummy.push_back(i+1);
	}
	dummy.push_back(n); // the dummy child
	scopes.push_back(dummy);

	// Create an undirected graph and check if disconnected.
	graph g(n+1);
	g.init(scopes);
	if (g.is_disconnected()) {
		return false;
	}

	// Create an undirected graph and check if disconnected.
	// graph g(n);
	// g.init(scopes);
	// if (g.is_disconnected()) {
	// 	return false;
	// }

	// Create the interval factors including the dummy one
	for (size_t i = 0; i < scopes.size(); ++i) {
		std::vector<size_t>& scope = scopes[i]; // 'i' is the child var
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		// assert(i == (int)scope.back()); // safety check
		f.set_child(scope.back());
		//f.fill_random_bayes();
		f.fill_random();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			// a = 0;
			// b = std::max(b, 0.6);
	
			if (i < scopes.size() - 1) {
				a = 0.0; //std::max(0.0, a - m_epsilon);
				b = std::max(b, 0.6);
				//b = b; //std::min(1.0, b + m_epsilon);
			} else {
				a = 0.0;
				b = 1.0;
			}

			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

bool generator::make_random_lcn2(){
	// Declare variables.
	size_t n = m_num_nodes;
	size_t c = m_num_nodes + m_num_extras;

	// Create the graph structure as follows:
	// - first variable in the ordering is root
	// - select pairs (x,y) => P(y|x)	
	std::vector<std::vector<size_t> > scopes;
	std::vector<size_t> scope(1);
	scope[0] = 0;
	scopes.push_back(scope);

	size_t count = 1;
	while (count < c) {

		// Pick the child and parent
		int child = randi2(n), parent;
		do {
			parent = randi2(n);
		} while (parent == child);
		std::vector<size_t> scope(2);
		scope[0] = parent;
		scope[1] = child;
		scopes.push_back(scope);

		++count;
	}

	// Create an undirected graph and check if disconnected.
	graph g(n);
	g.init(scopes);
	if (g.is_disconnected()) {
		return false;
	}

	// Create the interval factors including the dummy one
	for (size_t i = 0; i < scopes.size(); ++i) {
		std::vector<size_t>& scope = scopes[i];
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		f.set_child(scope.back());
		f.fill_random_bayes();
		//f.fill_random();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			if (k % 2 == 0) {
				a = std::max(0.0, a - m_epsilon);
				b = std::min(1.0, b + m_epsilon);
			} else {
				a = 0.0;
				b = 1.0;
			}

			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

bool generator::make_random_lcn3() {

	// Declare variables.
	size_t n = m_num_nodes;
	size_t c = m_num_nodes + m_num_extras;

	// Create a random ordering of the variables.
	std::vector<size_t> ordering(n);
	std::vector<size_t> position(n);
	for (size_t i = 0; i < n; ++i) {
		ordering[i] = i;
		position[i] = i;
	}

	// Randomly, switch pairs of variables.
	for (size_t i = 0 ; i < n ; ++i) {
		size_t j = randi2(n);
	
		// Switch variable i and j.
		int k = ordering[j];
		ordering[j] = ordering[i];
		ordering[i] = k;

		position[ordering[j]] = j;
		position[ordering[i]] = i;
	}


	// Create the graph structure as follows:
	// - select P(x_i+1|x_i)
	// - select additional random pairs (x,y) => P(y|x)	
	std::vector<std::vector<size_t> > scopes;
	std::vector<size_t> scope(1);
	scope[0] = ordering[0];
	scopes.push_back(scope);

	size_t count = 1;
	for (size_t i = 1; i < n; ++i) {
		int x = ordering[i - 1];
		int y = ordering[i];

		std::vector<size_t> scope(2);
		scope[0] = x; // parent
		scope[1] = y; // child
		scopes.push_back(scope);

		++count;
	}

	for (size_t j = 0; j < m_num_extras; ++j) {

		// Pick the child and parent
		int x = randi2(n), y;
		do {
			y = randi2(n);
		} while (y == x);
		std::vector<size_t> scope(2);
		scope[0] = ordering[x];
		scope[1] = ordering[y];
		scopes.push_back(scope);

		++count;
	}

	// Create an undirected graph and check if disconnected.
	graph g(n);
	g.init(scopes);
	if (g.is_disconnected()) {
		return false;
	}

	// Create the interval factors including the dummy one
	for (size_t i = 0; i < scopes.size(); ++i) {
		std::vector<size_t>& scope = scopes[i];
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		f.set_child(scope.back());
		f.fill_random_bayes();
		//f.fill_random();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			if (k % 2 == 0) {
				a = std::max(0.0, a - m_epsilon);
				b = std::min(1.0, b + m_epsilon);
			} else {
				a = 0.0;
				b = 1.0;
			}

			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

// Generate a grid network.
bool generator::make_grid() {
	
	// We assume square grids n-by-n = num_nodes
	size_t n_rows = (size_t) sqrt(m_num_nodes);
	size_t n_cols = n_rows;

	// Create the grid structure
	std::vector<std::vector<size_t> > scopes(m_num_nodes);
	for (size_t i = 0; i < n_rows; ++i) {
		size_t row_start = i * n_rows;

		// first row
		if (0 == i) {
			for (size_t j = 0; j < n_cols; ++j) {
				// first row, first column
				if (0 == j) {
					std::vector<size_t> scope(1);
					scope[0] = 0;

					scopes[0] = scope;
					continue;
				}

				// first row, 2-M'th column
				std::vector<size_t> scope(2);

				size_t c = row_start + j;
				scope[0] = c - 1;
				scope[1] = c;
				
				scopes[c] = scope;
			}

			continue;
		}

		// Rows 2-N'th
		for (size_t j = 0; j < n_cols; ++j) {
			// other row, first column
			if (0 == j) {
				std::vector<size_t> scope(2);
				
				size_t c = row_start;
				scope[0] = c - n_rows;
				scope[1] = c;

				scopes[c] = scope;				
				continue;
			}

			// Other row, columns 2-M'th
			std::vector<size_t> scope(3);

			size_t c = row_start + j ;
			scope[0] = c - 1;
			scope[1] = c - n_rows;
			scope[2] = c;

			scopes[c] = scope;
		}
	}
	
	// Create the interval factors 
	for (size_t i = 0; i < scopes.size(); ++i) {
		std::vector<size_t>& scope = scopes[i]; // 'i' is the child var
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		assert(i == (int)scope.back()); // safety check
		f.set_child(scope.back());
		f.fill_random_bayes();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			a = std::max(0.0, a - m_epsilon);
			b = std::min(1.0, b + m_epsilon);
			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

// Generate a ktree network
bool generator::make_ktree() {

	size_t n = m_num_nodes;
	size_t k = m_ksize;

	// A temporary variable for storing the graph
	std::vector<std::vector<size_t> > g;
	g.resize(n);
	for (size_t i = 0; i < n; i++) {
		g[i].resize(n);
	}

	// Initialize the graph to contain no edges
	// If there is an edge from vertex i to j then graph[i][j] =1
	// else graph[i][j]=0
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			g[i][j] = 0;
		}
	}

	std::vector<size_t> new_scope;
	std::vector<std::vector<size_t> > scopes(n);
	std::vector<std::vector<size_t> > cliques;
	// The way we generate k-tree is by first forming a clique of size k+1
	// Then taking a new vertex and making it adjacent to a clique of size k.

	// Form the first clique (also the root scopes)
	std::vector<size_t> new_clique;
	for (size_t i = 0; i < k; i++) {
		new_clique.push_back(i);
		new_scope.clear();
		new_scope.push_back(i);
		scopes[i] = new_scope;
	}
	cliques.push_back(new_clique);

	// Form the remaining cliques
	for (size_t i = k; i < n; i++) {
		size_t num_cliques = cliques.size();
		int cli = randi2(num_cliques);
		
		// Add a new scope with "i" as child and "cli" as parents.
		new_scope.clear();
		new_scope.push_back(i);
		copy(cliques[cli].begin(), cliques[cli].end(), back_inserter(new_scope));
		scopes[i] = new_scope;

		// Generate new k-cliques containing newly added "i"
		for (size_t j = 0; j < k; j++) {
			new_clique = std::vector<size_t>();
			for (size_t t = 0; t < k; t++)
				if (t == j) {
					new_clique.push_back(i);
				} else {
					new_clique.push_back(cliques[cli][t]);
				}
			cliques.push_back(new_clique);
		}
	}

	// Create the directed graph corresponding to the scopes.
	std::vector<std::vector<size_t> >::iterator si = scopes.begin();
	for (; si != scopes.end(); ++si) {
		std::vector<size_t>& sc = (*si);
		for (size_t ii = 1; ii < sc.size(); ++ii) {
			g[sc[ii]][sc[0]] = 1;	// add edge from "parent" to "child"
		}
	}

	typedef struct {
		int v1;
		int v2;
	} Edge;

	// Create a list of directed edges.
	int e = 0;
	std::vector<Edge> edges;
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++) {
			if (g[i][j] == 1) {
				e++;
				Edge edge;
				edge.v1=i;
				edge.v2=j;
				edges.push_back(edge);
			}
		}
	}

	// Remove a percent of the edges. The resulting graph
	// is a partial directed k-tree.
	int count = e*(m_kpercent)/100;
	e = e-count;
	for (size_t i = 0; i < count; i++) {
		bool dup;
		do {
			dup = false;
			int temp = randi2((long)edges.size());
			if (edges[temp].v1 < 0) {
				dup = true;
			} else {
				g[edges[temp].v1][edges[temp].v2] = 0;
				edges[temp].v1 = -1;
				edges[temp].v2 = -1;
			}
		}
		while(dup);
	}

	std::vector<std::vector<size_t> > cpts;
	for (si = scopes.begin(); si != scopes.end(); ++si)	{
		std::vector<size_t> &sc = (*si);
		std::vector<size_t> tmp; 
		tmp.push_back(sc[0]);
		for (size_t ii = 1; ii < sc.size(); ++ii)
			if (1 == g[sc[ii]][sc[0]])
				tmp.push_back(sc[ii]);

		size_t p = 0;
		std::vector<size_t> argv(tmp.size());
		for (size_t ii = 1; ii < tmp.size(); ++ii)
			argv[p++] = tmp[ii];
		argv[p] = tmp[0];

		cpts.push_back(argv);
	}

	// Create the interval factors 
	for (size_t i = 0; i < cpts.size(); ++i) {
		std::vector<size_t>& scope = cpts[i]; // 'i' is the child var
		variable_set argv;
		for (size_t j = 0; j < scope.size(); ++j) {
			argv |= variable(scope[j], 2);
		}

		interval f = interval(argv);
		assert(i == (int)scope.back()); // safety check
		f.set_child(scope.back());
		f.fill_random_bayes();
		for (size_t k = 0; k < f.numel(); ++k) {
			interval::value val = f[k];
			double a = val.first;
			double b = val.second;
			a = std::max(0.0, a - m_epsilon);
			b = std::min(1.0, b + m_epsilon);
			f[k] = interval::value(a, b);
		}

		m_factors.push_back(f);
	}

	fixup();

	return true;
}

/// Write the solution to the output stream
void generator::write_solution(std::ostream& os, int output_format) {

	assert(output_format == MERLIN_OUTPUT_UAI);
	std::cout << "[GEN] Writing converted credal net to the output file" << std::endl;
	
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
	// assert(n == m); // For Bayes nets and credal nets only
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

