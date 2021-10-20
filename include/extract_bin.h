#ifndef MGEMS_EXTRACT_BIN_H
#define MGEMS_EXTRACT_BIN_H

#include <vector>

#include "cxxio.hpp"

namespace mGEMS {
// mGEMS::ExtractBin
//   Extracts the reads assigned to a specific bin from the .fastq files.
//   Input:
//     `bin_assignments`: vector containing the ids of the reads assigned to this bin.
//     `in_strands`: the input .fastq files (e. g. forward and reverse strands).
//     `out_strands`: the output .fastq files.
void ExtractBin(const std::vector<uint32_t> &bin_assignments,
		std::vector<cxxio::In> &in_strands,std::vector<cxxio::Out> *out_strands);

// mGEMS::ReadBin
//   Reads in a bin that has been written to a file with mGEMS::WriteBin.
//   Input:
//     `stream`: Stream pointing to the bin file.
std::vector<uint32_t> ReadBin(std::istream &stream);
}

#endif
