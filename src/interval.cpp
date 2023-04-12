#include "interval.h"
#include "index.h"

namespace merlin {

interval::value interval::get_value(std::map<size_t, size_t>& config) {
    assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    return t_.at(i);
}

void interval::set_value(std::map<size_t, size_t>& config, value val) {
    assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    t_[i] = val;
}

potential interval::to_potential() {
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

        pot.add_p(f0);
        pot.add_q(f0);
        pot.add_p(f1);
        pot.add_q(f1);
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

            pot.add_p(f);
            pot.add_q(f);
		} // end while
    }

    return pot;
}

} // end namespace
