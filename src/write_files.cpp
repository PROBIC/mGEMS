#include <fstream>
#include <set>
#include <iostream>

#include "write_files.h"

void write_reads(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<std::string> &refnames,
		 const std::string &outfile) {

  std::cout << "n_refs: " << refnames.size() << std::endl;
  std::cout << "n_probs: " << probs.size() << std::endl;
  std::vector<std::ofstream> ofs(refnames.size());
  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].open(outfile + "/" + refnames[i] + "_reads.txt");
  }

  std::cout << "writing reads..." << std::endl;
  for (auto kv : reads_in_ec) {
    bool any_probs = (probs.find(kv.first) != probs.end());
    if (any_probs) {
      for (size_t i = 0; i < refnames.size(); ++i) {
	if (probs.at(kv.first)[i]) {
	  for (auto read : kv.second) {
	    ofs[i] << read << '\n';
	  }
	}
      }
    }
  }

  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].close();
  }
}
