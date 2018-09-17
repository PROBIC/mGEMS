#include <fstream>
#include <set>
#include <iostream>

#include "write_files.h"

void write_reads(const std::unordered_map<std::vector<bool>, long unsigned> &ec_to_id,
		 const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<std::string> &refnames,
		 const std::string &outfile) {

  std::cout << "n_refs: " << refnames.size() << std::endl;
  std::cout << "n_probs: " << probs.size() << std::endl;
  std::vector<std::ofstream> ofs(refnames.size());
  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].open(outfile + "/" + refnames[i] + "_reads.txt");
  }

  // std::ofstream control_file(outfile + "/" + "control.txt");
  // for (auto kv : reads_in_ec) {
  //   if (ec_to_id.find(kv.first) != ec_to_id.end()) {
  //     control_file << "reads in ec " << ec_to_id.at(kv.first) << ": " << kv.second.size() << '\n';
  //   }
  // }
  // control_file.close();

  std::cout << "writing reads..." << std::endl;
  for (auto kv : ec_to_id) {
    for (size_t i = 0; i < refnames.size(); ++i) {
      bool reads_exist = (reads_in_ec.find(kv.first) != reads_in_ec.end());
      bool probs_exist = (probs.find(kv.second) != probs.end());
      if (reads_exist && probs_exist) {
	if (probs.at(kv.second)[i] == true) {
	  for (size_t j = 0; j < reads_in_ec.at(kv.first).size(); ++j) {
	    ofs[i] << reads_in_ec.at(kv.first)[j] << '\n';
	  }
	}
      }
    }
  }

  for (size_t i = 0; i < refnames.size(); ++i) {
    ofs[i].close();
  }
}
