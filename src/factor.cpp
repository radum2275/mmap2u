#include "factor.h"
#include "index.h"

namespace merlin {

factor::value factor::get_value(std::map<size_t, size_t>& config) {
    assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    return t_.at(i);
}

void factor::set_value(std::map<size_t, size_t>& config, value val) {
    assert(config.size() == v_.size());
    config_index idx(v_, true); // default big endian
    size_t i = idx.convert(config);
    assert(i >= 0 && i < t_.size());
    t_[i] = val;
}

} // end namespace
