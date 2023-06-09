/*
 * algorithm.h
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

/// \file algorithm.h
/// \brief Algorithm interface
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#ifndef IBM_MERLIN_ALGORITHM_H_
#define IBM_MERLIN_ALGORITHM_H_

#include "base.h"
#include "interval.h"

///
/// \namespace merlin 
/// \brief The default namespace of the library.
///
namespace merlin {

///
/// \brief Interface for all inference algorithms.
///
/// The class provides the functional interface for any probabilistic
/// inference algorithm. It is a pure virtual class and therefore it should
/// be implemented by a concrete (derived) inference algorithm (e.g., wmb).
///
class algorithm {
public:

	///
	/// \brief Constructs an empty algorithm.
	///	
	algorithm() :
		m_stop_iter(0), m_stop_obj(0),
		m_stop_msg(0), m_start_time(0) {};

	///	
	/// \brief Destructs the algorithm.
	///
	virtual ~algorithm() {}; 
	
	///
	/// \brief Inialize the algorithm.
	///
	/// The initialization step must be performed before running the algorithm.
	///
	virtual void init() = 0;

	///
	/// \brief Run the algorithm.
	///
	/// Solves a probabilistic inference task.
	///
	virtual void run() = 0;

	// For variational methods (the concept of "belief" may vary):

	///
	/// \brief Belief associated with a node of the graph.
	///
	/// \fn const factor& belief(size_t i) const = 0
	/// \return the belief associated with a node of the graph. 
	///		Nodes correspond typically to variables. The belied is a function
	///		represented by the Factor class.
	/// \param i 	The index of the node in the graph (from 0)
	///
	virtual const interval& belief(size_t i) const = 0;

	///
	/// \brief Belief associated with a particular variable.
	///
	/// \param v 	The variable to compute the belief of
	/// \return the belief associated with a particular variable in the
	///		graphical model. The belief is a function
	///		represented by the Factor class.
	///
	virtual const interval& belief(variable v) const = 0;

	///
	/// \brief Beliefs associated with the variables.
	///
	/// \fn const std::vector<Factor>& beliefs() const = 0
	/// \return the beliefs associated with the variables in the
	///		graphical model (one belief for each of the variables). The output
	///		vector is indexed by the same indexes used for the variables.
	///
	virtual const std::vector<interval>& beliefs() const = 0;	// or all of them


	///
	/// \brief Set the number of iterations.
	///
	/// Stop when d*(# factors) updates have been done.
	///
	void set_stop_iter(double d) {
		m_stop_iter = d;
	}

	///
	/// \brief Set the minumum objective change before stopping.
	///
	/// Stop when objective change is less than *d*.
	/// \param d 	The objective change
	///
	void set_stop_obj(double d) {
		m_stop_obj = d;
	}

	///
	/// \brief Set the minimum ammount a message changes before stopping.
	///
	/// Stop when message updates are less than *d*.
	/// \param d 	The tolerated message update
	///
	void set_stop_msg(double d) {
		m_stop_msg = d;
	}

protected:
	// Properties such as # of iterations, cpu time, delta objective

	double m_stop_iter;					///< Stopping criteria (iterations)
	double m_stop_obj;					///< Stopping criteria (objective change)
	double m_stop_msg;	 				///< Stopping criteria (message change)
	double m_start_time;				///< Start time

};

} // namespace

#endif /* IBM_MERLIN_ALGORITHM_H_ */
