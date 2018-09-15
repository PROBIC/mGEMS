#include <fstream>
#include <set>

#include "write_files.h"

void write_reads(const std::map<std::vector<bool>, long unsigned> &ec_to_id, const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec, std::unordered_map<long unsigned, std::vector<bool>> probs, std::vector<std::string> refnames, const std::string &outfile) {

  std::vector<std::ofstream> ofs(refnames.size());
  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].open(outfile + "/" + refnames[i]);
  }

  for (auto kv : ec_to_id) {
    for (size_t i = 0; i < refnames.size(); ++i) {
      if (probs.at(kv.second)[i] == true) {
	for (size_t j = 0; j < reads_in_ec.at(kv.first).size(); ++j) {
	  ofs[i] << reads_in_ec.at(kv.first)[j];
	}
      }
    }
  }

  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].close();
  }
}
