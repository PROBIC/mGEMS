#include <fstream>
#include <set>
#include <iostream>
#include <memory>

#include "assign_reads/write_files.h"
#include "zstr/zstr.hpp"

void write_reads(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<std::string> &refnames,
		 const std::vector<short unsigned> &group_indices,
		 const std::string &outfile,
		 const bool gzip_output) {

  std::cout << "n_refs: " << refnames.size() << std::endl;
  std::cout << "n_probs: " << probs.size() << std::endl;
  std::vector<std::unique_ptr<std::ostream>> ofs(refnames.size());
  for (auto i : group_indices) {
    std::string fname = outfile + "/" + refnames[i] + "_reads.txt";
    if (gzip_output) {
      ofs[i] = std::unique_ptr<std::ostream>(new zstr::ofstream(fname + ".gz"));
    } else {
      ofs[i] = std::unique_ptr<std::ostream>(new std::ofstream(fname));
    }
  }

  std::cout << "writing reads..." << std::endl;
  for (auto kv : reads_in_ec) {
    bool any_probs = (probs.find(kv.first) != probs.end());
    if (any_probs) {
      for (auto i : group_indices) {
	if (probs.at(kv.first)[i]) {
	  for (auto read : kv.second) {
	    *ofs[i] << read << '\n';
	  }
	}
      }
    }
  }

  for (auto i : group_indices) {
    ofs[i]->flush();
  }
}

void write_ecs(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec, const std::string &outfile, const bool gzip_output) {
  std::unique_ptr<std::ostream> of;
  std::string outfile_name = outfile + "/" + "ec_to_read.csv";
  if (gzip_output) {
    of = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile_name + ".gz"));
  } else {
    of = std::unique_ptr<std::ostream>(new std::ofstream(outfile_name));
  }
  if (of->good()) {
    for (auto kv : reads_in_ec) {
      *of << kv.first;
      for (auto read : kv.second) {
	*of << ',' << read;
      }
      *of << '\n';
    }
  }
  of->flush();
}
