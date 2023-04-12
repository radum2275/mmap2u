/*
 * bn2cn.h
 *
 *  Created on: 19 Apr 2021
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

/// \file bn2cn.h
/// \brief Convert a Bayes net to an interval credal net
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_BN2CN_H_
#define IBM_LOOPY_BN2CN_H_

#include "bayes_net.h"
#include "credal_net.h"
#include "algorithm.h"

namespace merlin {


/**
 * BN2CN
 *
 * Convert a Bayes net into an interval credal net
 *
 */

class bn2cn: public credal_net, public algorithm {
public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	bn2cn() : credal_net() {
		set_properties();
	}

	///
	/// \brief Constructor with a credal net.
	///
	bn2cn(const credal_net& cn) : credal_net(cn) {
		set_properties();
	}

	///
	/// \brief Destructor
	///
	~bn2cn() {
	};

	inline const interval& belief(size_t i) const {
		throw std::runtime_error("not implemented");
	}
	inline const interval& belief(variable v) const {
		throw std::runtime_error("not implemented");
	}
	inline const std::vector<interval>& beliefs() const {
		throw std::runtime_error("not implemented");
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
	MER_ENUM( Property , Epsilon,Verbose,Seed );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Epsilon=0.0,Verbose=1,Seed=0");
			return;
		}
		m_verbose = 1;
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Epsilon:
				m_epsilon = atof(asgn[1].c_str());
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
	/// \brief Initialize the Loopy 2U algorithm.
	///
	void init();

protected:
	// Members:

	double m_epsilon;								///< Epsilon for intervals
	size_t m_verbose;								///< Verbosity level
	size_t m_seed;									///< Random number generator seed
};

} // namespace




#endif /* IBM_MERLIN_CTE_H_ */
