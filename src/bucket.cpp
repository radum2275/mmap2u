#include "bucket.h"

namespace merlin {

    std::vector<potential> bucket::create_partition(int ibound) {
    
        std::vector<potential> partition;
    
        // Sort the potentials in decreasing order of scope size
        std::sort(m_potentials.begin(), m_potentials.end(), compare_potentials_desc);
    
        // Create the partition
        for (size_t i = 0; i < m_potentials.size(); ++i) {
            potential& pot = m_potentials[i];
            variable_set vs = pot.vars();
            bool found = false;
    
            for (size_t j = 0; j < partition.size(); ++j) {
                potential& mb = partition[j];
                variable_set vs2 = mb.vars();
                variable_set tmp = vs | vs2;
                if (tmp.nvar() <= ibound) {
                    mb.multiply(pot);
                    found = true;
                    break; // found the mini-bucket
                }
            }
    
            // if no mini-bucket found, create a new one
            if (!found) {
                potential new_mb(1.0);
                new_mb.multiply(pot);
                partition.push_back(new_mb);
            }
        }
        
        return partition;
    }
std::vector<potential> bucket::create_partition(int ibound, size_t query_type, 
    size_t potential_approx, size_t potential_size, double eps) {

    std::vector<potential> partition;

    // Sort the potentials in decreasing order of scope size
    std::sort(m_potentials.begin(), m_potentials.end(), compare_potentials_desc);

	// Create the partition
    for (size_t i = 0; i < m_potentials.size(); ++i) {
        potential& pot = m_potentials[i];
        variable_set vs = pot.vars();
        bool found = false;

        for (size_t j = 0; j < partition.size(); ++j) {
            potential& mb = partition[j];
            variable_set vs2 = mb.vars();
            variable_set tmp = vs | vs2;
            if (tmp.nvar() <= ibound) {
                mb.multiply(pot);
                found = true;
                mb.set_original(false);                
                break; // found the mini-bucket
            }
        }

        // if no mini-bucket found, create a new one
        if (!found) {
            potential new_mb(1.0);
            new_mb.multiply(pot);
            new_mb.set_original(false);
            partition.push_back(new_mb);
        }
    }
	
    // Approximate the potentials if required
    for (size_t i = 0; i < partition.size(); ++i) {
        potential& mb = partition[i];
        mb.approximate(potential_approx, potential_size, eps);
        if (query_type == MERLIN_MAP_MAXIMAX) {
            mb.maximize();
        } else if (query_type == MERLIN_MAP_MAXIMIN) {
            mb.minimize();
        }
    }

    return partition;
}

} // end namespace
