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

} // end namespace
