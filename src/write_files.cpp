#include <fstream>
#include <set>
#include <iostream>

#include "write_files.h"

void write_reads(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<std::string> &refnames,
		 const std::vector<short unsigned> &group_indices,
		 const std::string &outfile) {

  std::cout << "n_refs: " << refnames.size() << std::endl;
  std::cout << "n_probs: " << probs.size() << std::endl;
  std::vector<std::ofstream> ofs(refnames.size());
  for (auto i : group_indices) {
    ofs[i].open(outfile + "/" + refnames[i] + "_reads.txt");
  }

  std::cout << "writing reads..." << std::endl;
  for (auto kv : reads_in_ec) {
    bool any_probs = (probs.find(kv.first) != probs.end());
    if (any_probs) {
      for (auto i : group_indices) {
	if (probs.at(kv.first)[i]) {
	  for (auto read : kv.second) {
	    ofs[i] << read << '\n';
	  }
	}
      }
    }
  }

  for (auto i : group_indices) {
    ofs[i].close();
  }
}

void write_ecs(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec, const std::string &outfile) {
  std::ofstream of(outfile + "/" + "ec_to_read.csv");
  if (of.is_open()) {
    for (auto kv : reads_in_ec) {
      of << kv.first;
      for (auto read : kv.second) {
	of << ',' << read;
      }
      of << '\n';
    }
  }
  of.close();
}
