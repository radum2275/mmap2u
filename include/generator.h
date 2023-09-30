/*
 * generator.h
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

/// \file generator.h
/// \brief Generates credal networks
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_LOOPY_GENERATOR_H_
#define IBM_LOOPY_GENERATOR_H_

#include "bayes_net.h"
#include "credal_net.h"
#include "algorithm.h"

namespace merlin {

class generator : public credal_net, public algorithm {

public:
	typedef credal_net::findex findex;        ///< Factor index
	typedef credal_net::vindex vindex;        ///< Variable index
	typedef credal_net::flist flist;          ///< Collection of factor indices

public:

	///
	/// \brief Default constructor.
	///
	generator() : credal_net() {
		set_properties();
	}

	///
	/// \brief Constructor with a credal net.
	///
	generator(const credal_net& cn) : credal_net(cn) {
		set_properties();
	}

	///
	/// \brief Destructor
	///
	~generator() {
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
	/// \brief Run the generator.
	///
	void run();
	
	///
	/// \brief Properties of the algorithm
	///
	MER_ENUM( Property , Epsilon,Graph,Seed,Instances,Nodes,Parents,KSize,KPercent,Query,Samples,Extras,Evid );


	// Setting properties (directly or through property string):

	///
	/// \brief Set the properties of the algorithm.
	/// \param opt 	The string containing comma separated property value pairs
	///
	virtual void set_properties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			set_properties("Epsilon=0.0,Graph=random,Seed=0,Instances=10,Nodes=100,Parents=2,KSize=10,KPercent=25,Query=10,Samples=10,Extras=1,Evid=0");
			return;
		}
		std::vector<std::string> strs = merlin::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = merlin::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Epsilon:
				m_epsilon = atof(asgn[1].c_str());
			case Property::Graph: // rand,grid,ktree
				m_graph_type = asgn[1].c_str();
				break;
            case Property::Instances:
                m_num_instances = atol(asgn[1].c_str());
                break;
            case Property::Nodes:
                m_num_nodes = atol(asgn[1].c_str());
                break;
            case Property::Parents:
                m_num_parents = atol(asgn[1].c_str());
                break;
            case Property::KSize:
                m_ksize = atol(asgn[1].c_str());
                break;
            case Property::KPercent:
                m_kpercent = atol(asgn[1].c_str());
                break;
			case Property::Query:
				m_num_query = atol(asgn[1].c_str());
				break;
			case Property::Samples:
				m_num_samples = atol(asgn[1].c_str());
				break;
			case Property::Extras:
				m_num_extras = atol(asgn[1].c_str());
				break;
			case Property::Evid:
				m_num_evid = atol(asgn[1].c_str());
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
	/// \brief Initialize the generator.
	///
	void init();

	/// @brief Set the input instance filename
	/// @param filename name of the file
	void set_input_filename(std::string filename) {
		m_input_filename = filename;
	}

protected:
    bool make_random();
    bool make_grid();
    bool make_ktree();
	
	bool make_random_lcn();
	bool make_random_lcn2();
	bool make_random_lcn3();

	///
	/// \brief Cleanup
	///
	void cleanup();

protected:
	// Members:

	double m_epsilon;					///< Epsilon for intervals
	size_t m_seed;						///< Random number generator seed
    std::string m_graph_type;           ///< Problem type: random, grid, ktree
    size_t m_num_nodes;                 ///< Problem size in terms of number of variables
    size_t m_num_parents;               ///< Number of parents per CPT
    size_t m_num_instances;             ///< Number of instances to generate for given size
    size_t m_ksize;                     ///< Clique size for k-trees
    size_t m_kpercent;                  ///< Percentage of egdes removed from ktree
	size_t m_num_query;					///< Number of MAP query variables
	size_t m_num_samples;				///< Number of sample queries to generate
	size_t m_num_extras;				///< Number of extra LCN statements to generate
	size_t m_num_evid;					///< Number of evidence variables in LCN

	std::string m_input_filename;		///< Input instance filename (if any)	
	size_t m_query_perc;				///< Percentage of query vars
};

} // end namespace

#endif
