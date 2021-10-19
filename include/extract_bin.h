#ifndef MGEMS_EXTRACT_BIN_H
#define MGEMS_EXTRACT_BIN_H

#include <vector>

#include "cxxio.hpp"

namespace mGEMS {
void ExtractBin(const std::vector<uint32_t> &bin_assignments, std::vector<cxxio::In> &in_strands, std::vector<cxxio::Out> *out_strands);
std::vector<uint32_t> ReadBin(std::istream &stream);
}

#endif
