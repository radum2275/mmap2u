#include "interval.h"
#include "index.h"

namespace merlin {

interval::value interval::get_value(std::map<size_t, size_t>& config) {
    // assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    return t_.at(i);
}

void interval::set_value(std::map<size_t, size_t>& config, value val) {
    // assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    t_[i] = val;
}

// Convert an interval to a potential. Assumes bi-valued variables!
potential interval::to_potential(bool has_q) {
    potential pot;

    // Enumerate configurations of parents
    size_t child;
    std::vector<size_t> parents;
    for (variable_set::const_iterator it = v_.begin(); it != v_.end(); ++it) {
        variable v = (*it);
        if (v.label() != (size_t)c_) {
            parents.push_back(v.label());
        } else {
            child = v.label();
        }
    }

    if (parents.size() == 0) { // no parents
        std::map<size_t, size_t> config;
        config[child] = 1;
        value val = get_value(config);
        factor f0(v_);
        f0[0] = 1.0 - val.first;
        f0[1] = val.first;
        factor f1(v_);
        f1[0] = 1.0 - val.second;
        f1[1] = val.second;

        // Create the p-component
        pot.add_p(f0);
        pot.add_p(f1);

        // Create the q-component (if needed)
        if (has_q) {
            pot.add_q(f0);
            pot.add_q(f1);
        }

    } else { // parents
        
        // Generate all parents configs
        size_t n = parents.size();
        std::vector<int> vals(n, 0);
	    vals[n - 1] = -1;
        int i;

        std::vector<std::map<size_t, size_t> > parents_configs;
		while (true) {

			// Enumerate "parent" variables.
			for (i = n - 1; i >= 0; --i) {
				if (vals[i] < 1) break;
				vals[i] = 0;
			}

			if (i < 0) break;	// done;
			++vals[i];

			// NOW: all "parents" have a specific value combination.
			std::map<size_t, size_t> config;
            for (size_t j = 0; j < parents.size(); ++j) {
                config[parents[j]] = vals[j];
			}

            parents_configs.push_back(config);
		} // end while

        // Enumerate the extension of the credal set
        size_t m = parents_configs.size();
        std::vector<int> ivals(m, 0);
        ivals[m - 1] = -1;
		while (true) {

			// Enumerate "parent" variables.
			for (i = m - 1; i >= 0; --i) {
				if (ivals[i] < 1) break;
				ivals[i] = 0;
			}

			if (i < 0) break;	// done;
			++ivals[i];

			// NOW: all "parents configs" have a specific value combination.
            // for each configuration, there are 2 extremes of the credal set
            // enumerate all extreme points of the local convex sets (intervals)
            factor f(v_);
            for (size_t j = 0; j < parents_configs.size(); ++j) {
                std::map<size_t, size_t> config = parents_configs[j];
                config[child] = 1;
                value bounds = get_value(config);
                double prob = (ivals[j] == 0) ? bounds.first : bounds.second;
                f.set_value(config, prob);
                config[child] = 0;
                f.set_value(config, 1.0 - prob);
			}

            // Create the p-component
            pot.add_p(f);

            // Create the q-component
            if (has_q) {
                pot.add_q(f);
            }
		} // end while
    }

    pot.set_original(true);
    return pot;
}

// Convert an interval to a potential. Assumes multi-valued variables!
potential interval::to_potential_multi() {
    potential pot;

    // Enumerate configurations of parents
    int child;
    std::vector<size_t> parents;
    std::vector<size_t> parents_size;
    size_t child_size;
    for (variable_set::const_iterator it = v_.begin(); it != v_.end(); ++it) {
        variable v = (*it);
        if (v.label() != (size_t)c_) {
            parents.push_back(v.label());
            parents_size.push_back(v.states());
        } else {
            child = v.label();
            child_size = v.states(); // domain size of the child
        }
    }    

    if (parents.size() == 0) { // no parents
        std::vector<std::vector<int> > child_combinations;
        int K = child_size - 1;
        int N = child_size;
        std::string bitmask(K, 1); // K leading 1's
        bitmask.resize(N, 0); // N-K trailing 0's
        
        // Generate the k-1 combinations of values for the child
        do {
            // Store integers and permute bitmask
            std::vector<int> comb;
            for (int i = 0; i < N; ++i) { // [0..N-1] integers
                if (bitmask[i]) {
                    comb.push_back(i);
                }
            }
            child_combinations.push_back(comb);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

        // For each k-1 values combination, generate all low/high permutation
        for (size_t c = 0; c < child_combinations.size(); ++c) {
            std::vector<int>& comb = child_combinations[c];
            std::vector<bool> used(N, false);
            size_t left_out;
            for (size_t j = 0; j < comb.size(); ++j) {
                used[comb[j]] = true;
            }
            for (size_t j = 0; j < used.size(); ++j) {
                if (!used[j]) {
                    left_out = j;
                    break;
                }
            } 

            // Generate low/high combinations for each k-1 subset
            size_t n = comb.size();
            std::vector<int> vals(n, 0);
            vals[n - 1] = -1;
            int i;
    
            std::vector<std::map<size_t, size_t> > low_high_configs;
            while (true) {
    
                // Enumerate "parent" variables.
                for (i = n - 1; i >= 0; --i) {
                    if (vals[i] < 1) break;
                    vals[i] = 0;
                }
    
                if (i < 0) break;	// done;
                ++vals[i];
    
                // NOW: all "parents" have a specific value combination.
                // 0 - low, 1 - high
                std::map<size_t, size_t> config;
                for (size_t jj = 0; jj < comb.size(); ++jj) {
                    config[comb[jj]] = vals[jj];
                }
    
                low_high_configs.push_back(config);
            } // end while
            
            // Generate the factors
            for (size_t j = 0; j < low_high_configs.size(); ++j) {
                std::vector<double> distr(child_size, 0.0);
                double sum_probs = 0.0;
                std::map<size_t, size_t> scope;
                std::map<size_t, size_t>& config = low_high_configs[j];
                std::map<size_t, size_t>::iterator mi = config.begin();
                for (; mi != config.end(); ++mi) {
                    size_t k = mi->first;
                    size_t lh = mi->second;

                    scope[child] = k; // scope is instantiated
                    value val = get_value(scope); // get the corresp. bounds
                    double prob = (lh == 0 ? val.first : val.second);
                    distr[k] = prob;
                    sum_probs += prob;
                }

                // Safety checks
                if (sum_probs <= 1.0) {
                    // assert(sum_probs <= 1.0); // should be smaller than 1.0
                    scope[child] = left_out;
                    value val = get_value(scope);
                    double left_prob = 1.0 - sum_probs;
                    if (left_prob >= val.first && left_prob <= val.second) {
                        distr[left_out] = left_prob; // good distribution
                        factor f(v_, distr);
                        f.set_child(child);
                        pot.add_p(f); // add it to the p-component of the potential
                    }
                }
            }
        }
    } else { // parents
        
        // Generate all parents configs
        size_t n = parents.size();
        std::vector<int> vals(n, 0);
	    vals[n - 1] = -1;
        int i;

        std::vector<std::map<size_t, size_t> > parents_configs;
		while (true) {

			// Enumerate "parent" variables.
			for (i = n - 1; i >= 0; --i) {
                int last = parents_size[i] - 1;
                if (vals[i] < last) break;
                vals[i] = 0;
			}

			if (i < 0) break;	// done;
			++vals[i];

			// NOW: all "parents" have a specific value combination.
			std::map<size_t, size_t> config;
            for (size_t j = 0; j < parents.size(); ++j) {
                config[parents[j]] = vals[j];
			}

            parents_configs.push_back(config);
		} // end while

        // Enumerate the extension of the credal set
        // NOW: all "parents configs" have a specific value combination.
        // For each configuration of the parents, enumerate the low/high
        // points for all combinations of the child's domain values
        for (size_t j = 0; j < parents_configs.size(); ++j) {
            std::map<size_t, size_t> config = parents_configs[j];

            // Generate all combinations of the child variable
            std::vector<std::vector<int> > child_combinations;
            int K = child_size - 1;
            int N = child_size;
            std::string bitmask(K, 1); // K leading 1's
            bitmask.resize(N, 0); // N-K trailing 0's
            
            // Generate the k-1 combinations of values for the child
            do {
                // Store integers and permute bitmask
                std::vector<int> comb;
                for (int k = 0; k < N; ++k) { // [0..N-1] integers
                    if (bitmask[k]) {
                        comb.push_back(k);
                    }
                }
                child_combinations.push_back(comb);
            } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    
            // For each k-1 values combination, generate all low/high permutation
            for (size_t c = 0; c < child_combinations.size(); ++c) {
                std::vector<int>& comb = child_combinations[c];
                std::vector<bool> used(N, false);
                size_t left_out;
                for (size_t j = 0; j < comb.size(); ++j) {
                    used[comb[j]] = true;
                }
                for (size_t j = 0; j < used.size(); ++j) {
                    if (!used[j]) {
                        left_out = j;
                        break;
                    }
                } 
    
                // Generate low/high combinations for each k-1 subset
                size_t nn = comb.size();
                std::vector<int> child_vals(nn, 0);
                child_vals[nn - 1] = -1;
                // int i;
        
                std::vector<std::map<size_t, size_t> > low_high_configs;
                while (true) {
        
                    // Enumerate "parent" variables.
                    for (i = nn - 1; i >= 0; --i) {
                        if (child_vals[i] < 1) break;
                        child_vals[i] = 0;
                    }
        
                    if (i < 0) break;	// done;
                    ++child_vals[i];
        
                    // NOW: all "parents" have a specific value combination.
                    // 0 - low, 1 - high
                    std::map<size_t, size_t> child_config;
                    for (size_t jj = 0; jj < comb.size(); ++jj) {
                        child_config[comb[jj]] = child_vals[jj];
                    }
        
                    low_high_configs.push_back(child_config);
                } // end while
                
                // Generate the factors
                for (size_t jj = 0; jj < low_high_configs.size(); ++jj) {
                    std::vector<double> distr(child_size, 0.0);
                    double sum_probs = 0.0;
                    std::map<size_t, size_t> scope = config; // parents config
                    std::map<size_t, size_t>& lh_config = low_high_configs[jj];
                    std::map<size_t, size_t>::iterator mi = lh_config.begin();
                    for (; mi != lh_config.end(); ++mi) {
                        size_t k = mi->first;
                        size_t lh = mi->second;
    
                        scope[child] = k; // full config of the scope
                        value val = get_value(scope);
                        double prob = (lh == 0 ? val.first : val.second);
                        distr[k] = prob;
                        sum_probs += prob;
                    }

                    // Safety checks
                    if (sum_probs <= 1.0) {
                        // assert(sum_probs <= 1.0); // should be smaller than 1.0
                        scope[child] = left_out;
                        value val = get_value(scope);
                        double left_prob = 1.0 - sum_probs;
                        if (left_prob >= val.first && left_prob <= val.second) {
                            distr[left_out] = left_prob; // good distribution
                            factor f(v_, distr);
                            f.set_child(child);
                            pot.add_p(f); // add it to the p-component of the potential
                        }
                    }
                }
            }
        
		} // end while
    }

    pot.set_original(true);
    return pot;
}

} // end namespace
