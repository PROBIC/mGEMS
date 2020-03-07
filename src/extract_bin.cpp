#include "extract_bin.h"

#include <fstream>
#include <algorithm>

namespace mGEMS {
  void ProcessStrand(const std::vector<uint32_t> &bin_assignments, std::istream &instrand, std::ostream &outstrand) {
  std::string line;
  uint32_t line_nr = 0;
  uint32_t bin_id = bin_assignments.size();
  while (std::getline(instrand, line)) {
    ++line_nr;
    uint32_t read_id = (line_nr - 1)/4 + 1;
    if (bin_assignments[bin_id] == read_id) {
      outstrand << line << '\n';
      for (uint8_t j = 0; j < 3; ++j) {
	std::getline(instrand, line);
	++line_nr;
	outstrand << line << '\n';
      }
      --bin_id;
    } else {
      for (uint8_t j = 0; j < 3; ++j) {
	std::getline(instrand, line);
	++line_nr;
      }
    }
  }
  outstrand.flush();
}

void ExtractBin(const std::vector<uint32_t> &bin_assignments, std::ostream* out_strand_1, std::ostream* out_strand_2, std::istream *in_strand_1, std::istream *in_strand_2) {
  ProcessStrand(bin_assignments, *in_strand_1, *out_strand_1);
  ProcessStrand(bin_assignments, *in_strand_2, *out_strand_2);
}

std::vector<uint32_t> ReadBin(std::istream &stream) {
  std::vector<uint32_t> out;
  std::string line;
  while (std::getline(stream, line)) {
    out.emplace_back(std::stoul(line));
  }
  return out;
}
}
