/*
 * credal_net.h
 *
 *  Created on: Oct 6, 2020
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

/// \file credal_net.h
/// \brief A credal network (imprecise Bayes net)
/// \author Radu Marinescu radu.marinescu@ie.ibm.com


#ifndef IBM_LOOPY_CREDAL_NET_H_
#define IBM_LOOPY_CREDAL_NET_H_

#include "enum.h"
#include "factor.h"
#include "graph.h"
#include "directed_graph.h"

namespace merlin {

///
/// \brief Credal net base class.
///
/// Simplest form of a *credal net*, namely just a collection of 
/// variables and factors (defined on subsets of variables).
/// Internally, factors and variables are mainly referenced by integer indices, with
///  0 <= f < nFactors() and 0 <= v < nvar().
///
class credal_net: public directed_graph {
public:

	// Useful typedefs:

	typedef size_t findex;			///< Factor index
	typedef size_t vindex;			///< Variable index
	typedef my_set<findex> flist; 	///< Collection of factor indices

	// Constructors, copy and assignment:

	///
	/// \brief Creates an empty credal net.
	///
	credal_net() :
		m_factors(), m_vadj(), m_dims() {
	};

	///
	/// \brief Creates a credal net by copying an existing one.
	/// \param gm 	An object of the same type
	///
	credal_net(const credal_net& cn) :
			directed_graph((directed_graph&) cn), m_factors(cn.m_factors),
			m_vadj(cn.m_vadj), m_dims(cn.m_dims) {
	};

	///
	/// \brief Assignment operator (deep copy).
	/// \param cn 	The credal net to be copied from
	/// \return a reference to the modified object containing the copied content.
	///
	credal_net& operator=(const credal_net& cn) {
		directed_graph::operator=((directed_graph&) cn); // copy graph elements over
		m_factors = cn.m_factors;				// copy list of factors
		m_vadj = cn.m_vadj;						// variable dependencies
		m_dims = cn.m_dims;						// and variable dimensions over
		return *this;
	};

	///
	/// \brief Clone the credal net.
	/// \return the pointer to the new object representing the copy of 
	/// the current model.
	///
	credal_net* clone() {
		credal_net* cn = new credal_net(*this);
		return cn;
	};

	///
	/// \brief Constructor from a list of factors.
	/// \param fs 	The list of factors
	///
	credal_net(const std::vector<factor>& fs) :
			directed_graph(), m_factors(fs), m_vadj(), m_dims() {
		fixup();
	};

	///
	/// \brief Creates a credal net from input iterators.
	/// \param first 	The iterator to beginning
	/// \param last 	The iterator to end
	///
	template<class InputIterator>
	credal_net(InputIterator first, InputIterator last) :
			m_factors(first, last), m_vadj(), m_dims() {
		fixup();
	}

	///
	/// \brief Destroys the credal net.
	///
	virtual ~credal_net() {
	};

	///
	/// \brief Read the credal net from a file in the UAI format.
	///
	///	For details on the UAI file format see the main documentation of
	/// the library. In this format, the factor tables are assumed to be
	/// represented using the "least significant bit", namely the last variable
	/// in the factor scope changes the fastest. Moreover, the internal
	/// representation of the factor assumes that the scope is ordered
	/// lexicographically.
	///
	/// \param file_name 	The full path to the file
	///
	void read(std::istream& is) {

		// Read the header
		size_t nvar, ncliques, csize, v, nval;
		std::string st;
		is >> st;
		int input_format = MERLIN_UNKNOWN;
		if ( st.compare("INTERVAL") == 0 ) {
			input_format = MERLIN_INPUT_INTERVAL;
		} else if ( st.compare("BAYES") == 0 ) {
			input_format = MERLIN_INPUT_BAYES;
		} else {
			std::string err_msg("Merlin only supports the BAYES net or INTERVAL credal net file formats.");
			throw std::runtime_error(err_msg);
		}

		// Read the number of variables and their domains
		is >> nvar;
		std::vector<size_t> dims(nvar);
		for (size_t i = 0; i < nvar; i++)
			is >> dims[i];

		// Read the number of factors and their scopes (scope is a variable_set)
		is >> ncliques;
		std::vector<std::vector<variable> > cliques(ncliques);
		std::vector<variable_set> sets(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> csize;
			cliques[i].reserve(csize);
			for (size_t j = 0; j < csize; j++) {
				is >> v;
				variable V(v, dims[v]);
				cliques[i].push_back(V);
				sets[i] |= V;
			}
		}

		// Read the factor tables (ensure conversion to ordered scopes)
		double lo, ub, pval;
		std::vector<factor> factors(ncliques);
		for (size_t i = 0; i < ncliques; i++) {
			is >> nval;
			assert(nval == sets[i].num_states());
			factors[i] = factor(sets[i]); // preallocate memory
			factors[i].set_child(cliques[i].back().label());
			
			convert_index ci(cliques[i], false, true); // convert from source order (littleEndian) to target order (bigEndian)
			for (size_t j = 0; j < nval; j++) {
				size_t k = ci.convert(j);	// get the index in the factor table
				if (input_format == MERLIN_INPUT_INTERVAL) {
					is >> lo; // read the lower bound value
					is >> ub; // read the upper bound value
				} else if (input_format == MERLIN_INPUT_BAYES) {
					is >> pval;
					lo = ub = pval;
				}

				factors[i][k] = factor::value(lo, ub); // save the factor value into the table
			}
		}

		m_factors = factors;
		fixup();
	}

	// Basic accessors:

	///
	/// \brief Return the number of variables in the model.
	///
	size_t nvar() const {
		return m_vadj.size();
	};

	///
	/// \brief Convert a variable index to a *variable* object.
	/// \param i 	The variable index to be converted
	/// \return the *variable* object corresponding to the index.
	///
	variable var(vindex i) const {
		return variable(i, m_dims[i]);
	};

	///
	/// \brief Return the number of factors in the model.
	///
	size_t num_factors() const {
		return m_factors.size();
	};

	///
	/// \brief Accessor for a factor (const).
	/// \param idx 	The index of the factor
	/// \return the factor corresponding to that index (const reference).
	///
	const factor& get_factor(findex idx) const {
		return m_factors[idx];
	};

	///
	/// \brief Acessor for a factor.
	/// \param idx The index of the factor
	/// \return the factor corresponding to the input index.
	///
	factor& get_factor(findex idx) {
		return m_factors[idx];
	}
	
	///
	/// \brief Set a factor
	/// \param i	The index of the factor
	/// \param f	The new factor
	void set_factor(findex i, const factor& f) {
		m_factors[i] = f;
	}

	///
	/// \brief Accessor for the factor container.
	/// \return the list of factors of the model.
	///
	const std::vector<factor>& get_factors() const {
		return m_factors;
	};

	// Basic variable-based queries:

	///
	/// \brief Factors depending on a variable.
	/// \param v 	The variable object
	/// \return the list of factor indexes containing the variable in their scopes.
	///
	const flist& with_variable(const variable& v) const {
		return m_vadj[_vindex(v)];
	}

	///
	/// \brief Union of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist with_variable_set(const variable_set& vs) const { //   or on all of a set of variables
		flist fs = with_variable(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			fs &= with_variable(vs[v]);
		return fs;
	}

	///
	/// \brief Intersection of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist intersects(const variable_set& vs) const {
		flist fs = with_variable(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			fs |= with_variable(vs[v]);
		return fs;
	}

	///
	/// \brief Factors depending on a variable.
	/// \param v 	The variable object
	/// \return the list of factor indexes containing the variable in their scopes.
	///	
	flist contains(const variable& v) const {
		return with_variable(v);
	}

	///
	/// \brief Union of factors depending on a set of variables.
	/// \param vs 	The set of variable object
	/// \return the list of factor indexes containing those variables in their scopes.
	///
	flist contains(const variable_set& vs) const {
		return with_variable_set(vs);
	}
	
	///
	/// \brief Factors contained by a set of variables.
	/// \param vs 	The set of variables
	/// \return the list of factor indexes contained by the set of variables.
	///
	flist contained_by(const variable_set& vs) const {
		flist fs2, fs = intersects(vs);
		for (size_t f = 0; f < fs.size(); f++)
			if (m_factors[fs[f]].vars() << vs)
				fs2 |= fs[f];
		return fs2;
	}

	///
	/// \brief Markov blanket of a variable.
	/// \param v 	The variable object
	/// \return the set of variables that form the Markov blanket 
	/// 	(ie, the variables that *v* may depend on).
	///
	variable_set markov_blanket(const variable& v) const {  // variables that v may depend on
		variable_set vs;
		const flist& nbrs = with_variable(v);
		for (flist::const_iterator f = nbrs.begin(); f != nbrs.end(); ++f)
			vs |= get_factor(*f).vars();
		vs /= v;
		return vs;
	}

	///
	/// \brief Markov blanket of a set of variables.
	/// \param vs 	The set of variables
	/// \return the union of Markov blankets corresponding to the 
	/// 	variables in the input set.
	///
	variable_set markov_blanket(const variable_set& vs) const {
		variable_set ret = markov_blanket(vs[0]);
		for (size_t v = 1; v < vs.size(); v++)
			ret |= markov_blanket(vs[v]);
		return ret;
	}

	///
	/// \brief Full adjacency matrix.
	/// \return the adjacency list associated with each variable in the model.
	///
	std::vector<variable_set> mrf() const {		// full variable-to-var adjacency
		std::vector<variable_set> vvs;
		for (size_t v = 0; v < nvar(); ++v)
			vvs.push_back(markov_blanket(var(v)));
		return vvs;
	}

	// Factor ("node") manipulation operations:

	///
	/// \brief Add a new factor to the model.
	/// \param F 	The factor to be added
	/// \return the index associated with the newly added factor.
	///
	virtual findex add_factor(const factor& F) {         // add a factor to our collection
		const variable_set& v = F.vars();
		findex use = add_node();
		//if (use>=nFactors()) _factors.push_back(F); else _factors[use]=F;
		if (use >= num_factors()) {
			if (m_factors.capacity() > num_factors())
				m_factors.push_back(F);
			else {                           // if we'd need to copy, do it manually
				std::vector<factor> tmp;
				tmp.reserve(2 * num_factors());
				tmp.resize(num_factors() + 1);
				for (size_t i = 0; i < m_factors.size(); ++i)
					tmp[i].swap(m_factors[i]);
				tmp[num_factors()] = F;
				m_factors.swap(tmp);
			}
		} else
			m_factors[use] = F;

		insert(m_vadj, use, v);
		if (m_dims.size() < m_vadj.size())
			m_dims.resize(m_vadj.size(), 0);
		for (variable_set::const_iterator i = v.begin(); i != v.end(); ++i) { // look up dimensions if required
			if (m_dims[_vindex(*i)] == 0)
				m_dims[_vindex(*i)] = i->states();// add if we haven't seen this var
			else if (m_dims[_vindex(*i)] != i->states())	//   or check it against our current states
				throw std::runtime_error(
						"Incompatible state dimension in added factor");
		}
		return use;                                  // return the factor index used
	}

	///
	/// \brief Remove a factor from the model.
	/// \param idx 	The index of the factor to be removed
	/// 
	virtual void remove_factor(findex idx) {        // remove a factor from the collection
		erase(m_vadj, idx, get_factor(idx).vars());		// remove from variable lists
		m_factors[idx] = factor();                       // empty its position
		remove_node(idx);								// and remove the node
	}

	///
	/// \brief Remove all factors.
	///
	void clear_factors() {
		m_factors.clear();
		m_vadj.clear();
		m_vadj.resize(nvar());
		directed_graph::clear();
	};

	///
	/// \brief Find the factor with the smallest scope.
	/// \param fl 	The reference to a list of factors
	/// \return the index of the smallest factor (scope-wise).
	/// 
	findex smallest(const flist& fl) {
		assert(fl.size() > 0);
		findex ret = *fl.begin();
		for (flist::const_iterator f = fl.begin(); f != fl.end(); ++f)
			if (get_factor(ret).nvar() > get_factor(*f).nvar())
				ret = *f;
		return ret;
	}

	///
	/// \brief Find the factor with the largest scope.
	/// \param fl 	The reference to a list of factors
	/// \return the index of the largest factor (scope-wise).
	///
	findex largest(const flist& fl) {
		assert(fl.size() > 0);
		findex ret = *fl.begin();
		for (flist::const_iterator f = fl.begin(); f != fl.end(); ++f)
			if (get_factor(ret).nvar() < get_factor(*f).nvar())
				ret = *f;
		return ret;
	}

	// Check graphical model properties:

	///
	/// \brief Check if a binary model.
	/// \return *true* if all variables have at most two values in 
	/// 	their domains. Otherwise return *false*.
	///
	bool is_binary() const {
		for (size_t i = 0; i < m_dims.size(); ++i) {
			if (m_dims[i] > 2)
				return false;
		}
		return true;
	}
	
	///
	/// \brief Check if a pairwise model.
	/// \return *true* if all factors involve at most two variables.
	/// 	Otherwise return *false*.
	///
	bool is_pairwise() const {
		for (size_t i = 0; i < num_factors(); ++i) {
			if (m_factors[i].nvar() > 2)
				return false;
		}
		return true;
	}

	// Distribution-based operators:


	// Ordering: variable (elimination) orders and factor orders

	///
	/// \brief Find the induced width (complexity) of an elimination order.
	/// \param order 	The variable elimination order
	/// \return the induced width of the elimination order.
	///
	size_t induced_width(const variable_order_t& order) const {

		size_t width = 0;
		std::vector<variable_set> adj = mrf();
		size_t n = order.size();
		graph g(n); // create the undirected graph (moral graph)
		for (size_t i = 0; i < adj.size(); ++i) {
			const variable_set& vi = adj[i];
			for (variable_set::const_iterator cj = vi.begin();
					cj != vi.end(); ++cj) {
				size_t j = _vindex(*cj);
				g.add_edge(i, j);
			}
		}

		std::vector<size_t> position(n);
		for (size_t i = 0; i < n; ++i) {
			position[order[i]] = i;
		}

		// eliminate variables and create induced edges
		for (size_t i = 0; i < n; ++i) {
			size_t var = order[i];
			size_t pos = position[var];

			// find the neighbors appearing later in the ordering
			const my_set<edge_id> ns = g.neighbors(var);
			std::set<size_t> S;
			for (my_set<edge_id>::const_iterator si = ns.begin();
					si != ns.end(); ++si) {
				size_t j = si->second;
				if (position[j] > pos) {
					S.insert(j);
				}
			}

			// connect the neighbors appearing later in the ordering
			width = std::max(width, S.size());
			std::set<size_t>::iterator it1, it2;
			for (it1 = S.begin(); it1 != S.end(); ++it1) {
				it2 = it1;
				while (++it2 != S.end()) {
					size_t u = *it1, v = *it2;
					if (u != v) g.add_edge(u, v);
				}
			}
		}

		return width;
	}

    ///
    /// \brief Find a minfill variable elimination order.
    /// \return the variable ordering corresponding to the method, such that
    ///		the first variable in the ordering is eliminated first.
    ///
	variable_order_t order() const {
		variable_order_t order;
		order.resize(nvar());

		std::vector<variable_set> adj = mrf();
		typedef std::pair<double, size_t> node_score_t;
		typedef std::multimap<double, size_t> score_map_t;
		score_map_t scores;
		std::vector<score_map_t::iterator> reverse(nvar());

		for (size_t v = 0; v < nvar(); v++) { // get initial scores
			reverse[v] = scores.insert(node_score_t(order_score(adj, v), v));
		}

		for (size_t ii = 0; ii < nvar(); ++ii) { // iterate through, selecting variables
			score_map_t::iterator first = scores.begin(); // choose a random entry from among the smallest
			score_map_t::iterator last = scores.upper_bound(first->first);
			std::advance(first, randi(std::distance(first, last)));
			size_t i = first->second;

			order[ii] = var(i).label();  // save its label in the ordering
			scores.erase(reverse[i]);	 // remove it from our list
			variable_set vi = adj[i]; // go through adjacent variables (copy: adj may change)
			variable_set fix;		  //  and keep track of which need updating
			for (variable_set::const_iterator j = vi.begin(); j != vi.end(); ++j) {
				size_t v = _vindex(*j);
				adj[v] |= vi;             // and update their adjacency structures
				adj[v] /= var(i);
				if (fix.size() < scores.size()) {
					fix |= adj[v];	// come back and recalculate their scores
				}
			}
			for (variable_set::const_iterator j = fix.begin(); j != fix.end(); ++j) {
				size_t jj = j->label();
				scores.erase(reverse[jj]);	// remove and update (score,index) pairs
				reverse[jj] = scores.insert(node_score_t(order_score(adj, jj), jj));
			}
		}

		return order;
	}

    ///
    /// \brief Find a minfill variable elimination order.
    /// \param ord_type 	The ordering method
    /// \return the variable ordering corresponding to the method, such that
    ///		the first variable in the ordering is eliminated first.
    ///
	variable_order_t order2() const {

		// variable order to be computed
		variable_order_t order;

		// create a temporary graph
		graph G(nvar());
		std::vector<variable_set> adj = mrf();
		for (size_t i = 0; i < nvar(); ++i) {
			variable_set& vs = adj[i];
			for (variable_set::const_iterator vi = vs.begin();
					vi != vs.end(); ++vi) {
				size_t j = vi->label();
				if (j != i) G.add_edge(i, j);
			}
		}

		// initialize the order
		size_t n = nvar();
		order.reserve(n);

		// keeps track of node scores
		std::vector<int> scores(n);
		for (size_t i = 0; i < n; ++i) {
			scores[i] = order_score(G, i);
		}

		// eliminate nodes until all gone
		int min_score = -1, num_nodes = n;
		size_t width = 0;
		while (num_nodes != 0) {

			// keeps track of minimal score nodes
			std::vector<size_t> candidates; // minimal score of 1 or higher
			std::vector<size_t> simplicial; // simplicial nodes (score 0)

			min_score = std::numeric_limits<int>::max();

			// find node to eliminate
			for (size_t i = 0; i < n; ++i) {
				if (scores[i] == 0) { // score 0
					simplicial.push_back(i);
				} else if (scores[i] < min_score) { // new, lower score (but greater 0)
					min_score = scores[i];
					candidates.clear();
					candidates.push_back(i);
				} else if (scores[i] == min_score) { // current min. found again
					candidates.push_back(i);
				}
			}

			// eliminate all nodes with score=0 -> no edges will have to be added
			for (std::vector<size_t>::iterator it = simplicial.begin();
					it != simplicial.end(); ++it) {
				size_t v = (*it);
				order.push_back(v);
				--num_nodes;
				my_set<edge_id> temp = G.neighbors(v);
				width = std::max(width, temp.size());
				for (my_set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					size_t i = ci->first, j = ci->second;
					G.remove_edge(i, j);
				}
				G.remove_node(v); // and adj edges
				scores[v] = std::numeric_limits<int>::max();
			}

			// anything left to eliminate? If not, we are done!
			if (min_score == std::numeric_limits<int>::max()) {
				break;
			}

			// Pick one of the minimal score nodes (with score >= 1),
			// breaking ties randomly
			size_t cand = candidates[randi(candidates.size())];
			//size_t cand = candidates[0]; // first candidate in list (lexicograhically)
			order.push_back(cand);
			--num_nodes;

			// remember it's neighbors, to be used later
			my_set<edge_id> nlist = G.neighbors(cand);
			std::set<size_t> neighbors;
			for (my_set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				if (ci->first == cand && ci->second != cand)
					neighbors.insert(ci->second);
			}

			// connect neighbors in primal graph
			std::set<size_t>::iterator it1, it2;
			for (it1 = neighbors.begin(); it1 != neighbors.end(); ++it1) {
				it2 = it1;
				while (++it2 != neighbors.end()) {
					size_t i = *it1, j = *it2;
					G.add_edge(i, j);
				}
			}

			// compute candidates for score update (node's neighbors and their neighbors)
			width = std::max(width, neighbors.size());
			std::set<size_t> to_fix(neighbors);
			for (std::set<size_t>::const_iterator it = neighbors.begin();
					it != neighbors.end(); ++it) {
				size_t v = (*it);
				const my_set<edge_id>& temp = G.neighbors(v);
				for (my_set<edge_id>::const_iterator ci = temp.begin();
						ci != temp.end(); ++ci) {
					if (ci->first == v && v != ci->second)
						to_fix.insert(ci->second);
				}
			}
			to_fix.erase(cand);

			// remove node from primal graph
			for (my_set<edge_id>::const_iterator ci = nlist.begin();
					ci != nlist.end(); ++ci) {
				size_t i = ci->first, j = ci->second;
				G.remove_edge(i, j);
			}
			G.remove_node(cand); // and adj edges
			scores[cand] = std::numeric_limits<int>::max(); // tag score

			// update scores in primal graph (candidate nodes computed earlier)
			for (std::set<size_t>::const_iterator it = to_fix.begin();
					it != to_fix.end(); ++it) {
				size_t v = (*it);
				scores[v] = order_score(G, v);
			}
		}

		return order;
	}

	// Helpful manipulation functions for other data:

	///
	/// \brief Add the scope of a factor to the adjacency list
	/// \param adj 	The adjacency list to be modified
	/// \param idx 	The index of the factor
	/// \param vs 	The scope of the factor
	///
	void insert(std::vector<flist>& adj, findex idx, const variable_set& vs) {
		if (vs.nvar() > 0 && adj.size() <= vs.rbegin()->label()) // if we need to, expand our set of variables
			adj.resize(vs.rbegin()->label() + 1); //   to be large enough to be indexed by label
		for (size_t i = 0; i < vs.nvar(); ++i)
			adj[vs[i]] |= idx;			//   and add factor to adj list
	}

	///
	/// \brief Remove the scope of a factor from the adjacency list
	/// \param adj 	The adjacency list to be modified
	/// \param idx 	The index of the factor
	/// \param vs 	The scope of the factor
	///
	void erase(std::vector<flist>& adj, findex idx, const variable_set& vs) {
		for (size_t i = 0; i < vs.nvar(); ++i)
			adj[vs[i]] /= idx; 		//   remove a factor from each var's adj list
	}

	// Simple optimum selection routines:

	///
	/// \brief Output operator.
	/// \param out 	The reference of an output stream
	/// \param gm 	The reference of a graphical model
	/// \return the reference of the modified output stream containing the 
	/// 	formatted content of the graphical model.
	///
	friend std::ostream& operator<<(std::ostream& out, const credal_net& gm) {
		out << "[DEBUG] Dumping credal net content with " << gm.nvar() << " variables and "
				<< gm.num_factors() << " factors: " << std::endl;
		for (size_t j = 0; j < gm.get_factors().size(); j++) {
			out << " " << j << " " << gm.get_factors()[j] << std::endl;
		}
		out << "[DEBUG] Dumping directed graph" << std::endl;
		out << " - num nodes: " << gm.num_nodes() << std::endl;
		out << " - num edges: " << gm.num_edges() << std::endl;
		for (size_t i = 0; i < gm.edges().size(); ++i) {
			out << gm.edges()[i] << std::endl;
		}
		return out;
	}

protected:

	///
	/// \brief Look up the index of a variable.
	/// \param v 	The variable object
	/// \return the index corresponding to the variable in the model.
	size_t _vindex(const variable& v) const {
		return v.label();
	}

	///
	/// \brief Mutable accessor of adjacency.
	/// \param v 	The variable object
	/// \return the list of factors containing the variable in their scope.
	///
	flist& _with_variable(const variable& v) {
		return m_vadj[_vindex(v)];
	}

	///
	/// \brief Internal function for ensuring consistency of the model.
	///
	void fixup() {
		size_t nVar = 0;
		m_vadj.clear();
		m_dims.clear();
		for (std::vector<merlin::factor>::iterator f = m_factors.begin();
				f != m_factors.end(); ++f) {
			if (f->nvar())
				nVar = std::max(nVar, f->vars().rbegin()->label() + 1);
		}

		add_nodes(nVar); // add the variable nodes in the directed graph
		m_vadj.resize(nVar); // reserve the variable adjaceny lists (factors)
		m_dims.resize(nVar); // make space for variable inclusion mapping
		for (size_t f = 0; f < m_factors.size(); ++f) {	// for each factor,
			const variable_set& v = m_factors[f].vars(); // save the variables' dimensions and
			int child = m_factors[f].get_child();
			assert(child >= 0);
			for (variable_set::const_iterator i = v.begin(); i != v.end(); ++i) { // index this factor as including them
				m_dims[_vindex(*i)] = i->states(); // check against current values???
				_with_variable(*i) |= f;

				if (*i != (size_t)child) {
					add_edge(*i, child);
				}
			}
		}
	}

	// Internal helper functions:

	///
	/// \brief Compute the score of a variable for the elimination order.
	/// \param adj 		The adjacency lists
	/// \param i 		The index of the variable
	/// \return the score of the variable.
	///
	double order_score(const std::vector<variable_set>& adj, size_t i) const {
		double s = 0.0;
		for (variable_set::const_iterator j = adj[i].begin(); j != adj[i].end(); ++j) {
			s += (adj[i] - adj[_vindex(*j)]).size();
		}

		return s;
	}

	int order_score(const graph& g, size_t v) const {
		int s = 0.0;
		const my_set<edge_id>& nlist = g.neighbors(v);
		std::set<size_t> S;
		for (my_set<edge_id>::const_iterator ci = nlist.begin(); ci != nlist.end(); ++ci) {
			size_t i = ci->first, j = ci->second;
			if (i == v && i != j)
				S.insert(j);
		}
		std::set<size_t>::iterator it1, it2;
		for (it1 = S.begin(); it1 != S.end(); ++it1) {
			it2 = it1;
			while (++it2 != S.end()) {
				size_t i = *it1, j = *it2;
				if (g.edge(i, j) == edge_id::NO_EDGE)
					++s;
			}
		}

		return s;
	}

	///
	/// \brief Create a random elimination order.
	///
	variable_order_t order_random() const {
		variable_order_t order;
		order.resize(nvar());
		for (size_t i = 0; i < nvar(); i++)
			order[i] = var(i).label();		// build a list of all the variables
		std::random_shuffle(order.begin(), order.end());// and randomly permute them
		return order;
	}

protected:
	// Members:

	std::vector<factor> m_factors;  ///< Collection of all factors in the model.

protected:
	// Members:

	std::vector<flist> m_vadj;		///< Variable adjacency lists (variables to factors)
	std::vector<double> m_dims;		///< Dimensions of variables as stored in graphical model object
};

} // namespace


#endif /* IBM_MERLIN_GRAPHICAL_MODEL_H_ */
