#ifndef MGEMS_EXTRACT_BIN_H
#define MGEMS_EXTRACT_BIN_H

#include <vector>

#include "file.hpp"

namespace mGEMS {
void ExtractBin(const std::vector<uint32_t> &bin_assignments, std::vector<File::In> &in_strands, std::vector<File::Out> *out_strands);
std::vector<uint32_t> ReadBin(std::istream &stream);
}

#endif
