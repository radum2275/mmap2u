/*
 * loopy2u.h
 *
 *  Created on: 06 Oct 2020
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

/// \file loopy2u.h
/// \brief Loopy2U for credal nets with intervals and binary variables
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_LOOPY2U_H_
#define IBM_LOOPY_LOOPY2U_H_

#include "credal_net.h"
#include "algorithm.h"

namespace merlin {


/**
 * Loopy 2U (L2U)
 *
 * Tasks supported: MAR
 *
 */

class loopy2u: public credal_net, public algorithm {
public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

	class message { // directed edge message Ui -> X
	public:
		variable parent; // parent (Ui)
		variable child; // child (X)
		interval pi; // message from parent to child: pi(Ui)
		interval lambda; // message from child to parent: scalar
		bool evidence; // dummy message for evidence nodes
		message(variable p, variable c) {
			parent = p;
			child = c;
			pi = interval(interval::value(1.0, 1.0)); // scalar (only p=1)
			lambda = interval(interval::value(1.0, 1.0)); // scalar
			evidence = false;
		}
	};

public:

	///
	/// \brief Default constructor.
	///
	loopy2u() : credal_net() {
		set_properties();
	}

	///
	/// \brief Constructor with a credal net.
	///
	loopy2u(const credal_net& cn) : credal_net(cn) {
		set_properties();
	}

	///
	/// \brief Constructor with a set of factors.
 	///
	loopy2u(const std::vector<interval>& fs) : credal_net(fs) {
		set_properties();
	}

	///
	/// \brief Destructor
	///
	~loopy2u() {
	};

	inline const interval& belief(size_t i) const {
		return m_beliefs[i];
	}
	inline const interval& belief(variable v) const {
		return m_beliefs[v];
	}
	inline const std::vector<interval>& beliefs() const {
		return m_beliefs;
	}

	///
	/// \brief Write the solution to the output stream.
	/// \param out		 		The output stream
	/// \param output_format	The output format (json or uai)
	///
	void write_solution(std::ostream& out, int output_format);

	///
	/// \brief Run the Loopy 2U algorithm.
	///
	void run();

	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , StopIter,Threshold,Verbose,Seed );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the variable order.
	///
	inline void set_order(const variable_order_t& ord) {
		m_order = ord;
	}

	///
	/// \brief Get the variable order.
	///
	inline const variable_order_t& get_order() const {
		return m_order;
	}

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("StopIter=10,Threshold=1e-06,Verbose=1,Seed=0");
			return;
		}
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::StopIter:
				m_iterations = atol(asgn[1].c_str());
				break;
			case Property::Threshold:
				m_threshold = atof(asgn[1].c_str());
				break;
			case Property::Verbose:
				m_verbose = atol(asgn[1].c_str());
				break;
			case Property::Seed:
				m_seed = atol(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	///
	/// \brief Reset the internal state of the solver
	///
	void reset();

	///
	/// \brief Set the evidence variables
	///
	void set_evidence(const std::map<size_t, size_t>& ev) {
		m_evidence = ev;
	}

	///
	/// \brief Initialize the Loopy 2U algorithm.
	///
	void init();

	///
	/// \brief Update the beliefs for all variables
	///
	void update_beliefs();

	///
	/// \brief Update the belief of a variable
	///
	void upate_belief(variable x);

	///
	/// \brief Compute the pi(x) message
	///
	/// \param x it the current variable
	///
	void pi(variable x);
	
	///
	/// \brief Compute the lambda(x) message
	///
	/// \param x it the current variable
	///	
	void lambda(variable x);
	
	///
	/// \brief Compute the pi(x-child) message
	///
	/// \param x it the current variable
	/// \param y is the child variable
	///
	void pi(variable x, variable y, message& m);

	///
	/// \brief Compute the lambda(x-parent) message
	///
	/// \param x it the current variable
	/// \param u is the parent variable
	///
	void lambda(variable x, variable u, message& m);

private:
	double hi(variable x, variable u, size_t ui, bool low, 
		std::vector<variable>& parents, std::vector<int>& extremes);
	double gi1(variable x, variable u, double L, 
		std::vector<variable>& parents, std::vector<int>& extremes);
	double gi2(variable x, variable u, double L, 
		std::vector<variable>& parents, std::vector<int>& extremes);

	void test();

protected:
	// Members:

	variable_order_t m_order;						///< Variable order
	std::vector<interval> m_beliefs; 					///< Marginals
	std::map<size_t, size_t> m_evidence;			///< Evidence
	size_t m_iterations;							///< Number of iterations
	double m_threshold;								///< Convergence threshold (default=1e-06)
	std::vector<size_t> m_schedule;					///< Propagation schedule
	size_t m_verbose;								///< Verbosity level
	std::vector<interval> m_pi;						///< Pi's for each variable x 
	std::vector<interval> m_lambda;					///< Lambda's for each variable x
	std::vector<loopy2u::message> m_messages;		///< Edge messages	
	std::vector<std::vector<size_t> > m_incoming;	///< Parents of a node (incoming)
	std::vector<std::vector<size_t> > m_outgoing;	///< Children of a node (outgoing)
	interval::value m_delta;							///< Average change in messages
	size_t m_seed;									///< Random number generator seed
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
