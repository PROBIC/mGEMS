#ifndef MGEMS_EXTRACT_BIN_H
#define MGEMS_EXTRACT_BIN_H

#include <vector>

#include "file.hpp"

namespace mGEMS {
void ExtractBin(const std::vector<uint32_t> &bin_assignments, std::ostream* out_strand_1, std::ostream* out_strand_2, std::istream* in_strand_1, std::istream* in_strand_2);
std::vector<uint32_t> ReadBin(std::istream &stream);
}

#endif
