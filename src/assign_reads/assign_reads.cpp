#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <cmath>

#include "arguments.h"
#include "assign_reads/read_files.h"
#include "assign_reads/write_files.h"
#include "zstr/zstr.hpp"

void multiply_abundances(std::vector<long double> &abundances, long double log_thresh) {
  for (size_t i = 0; i < abundances.size(); ++i) {
    abundances[i] = std::exp(std::log(abundances.at(i)) + log_thresh);
  }
}

int main(int argc, char* argv[]) {
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  std::cout << "Reading read assignments to equivalence classes" << std::endl;
  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::string assignments_path = std::string(GetOpt(argv, argv+argc, "-f"));
    zstr::ifstream assignments_file(assignments_path);
    read_assignments(assignments_file, &reads_to_ec);
  }

  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");

  std::string outfile_name = std::string(GetOpt(argv, argv+argc, "-o"));
  std::cout << "Reading abundances" << std::endl;
  std::string abundances_path = std::string(GetOpt(argv, argv+argc, "-a"));
  zstr::ifstream abundances_file(abundances_path);
  std::vector<std::string> ref_names;
  std::vector<long double> abundances = read_abundances(abundances_file, ref_names);
  double log_thresh = std::log1pl(-(long double)abundances.size()/(long double)reads_to_ec.size());
  if (CmdOptionPresent(argv, argv+argc, "-q")) {
    log_thresh += std::stold(std::string(GetOpt(argv, argv+argc, "-q")));
  }
  std::cout << std::exp(log_thresh) << std::endl;
  multiply_abundances(abundances, log_thresh);
    
  std::cout << "Reading probs" << std::endl;
  std::unordered_map<long unsigned, std::vector<bool>> probs;
  if (read_from_cin) {
    read_probs(abundances, std::cin, &probs);
  } else {
    std::string probs_path = std::string(GetOpt(argv, argv+argc, "-p"));
    zstr::ifstream probs_file(probs_path);
    read_probs(abundances, probs_file, &probs);
  }

  bool all_groups = !CmdOptionPresent(argv, argv+argc, "--groups");

  std::cout << "Assigning reads to reference groups" << std::endl;
  std::vector<short unsigned> group_indices;
  if (all_groups) {
    group_indices.resize(ref_names.size());
    for (size_t i = 0; i < ref_names.size(); ++i) {
      group_indices[i] = i;
    }
  } else {
    std::string groups_path = std::string(GetOpt(argv, argv+argc, "--groups"));
    std::ifstream groups_file(groups_path);
    read_groups(ref_names, groups_file, &group_indices);
  }
  std::vector<std::unique_ptr<std::ostream>> outfiles(ref_names.size());
  for (auto i : group_indices) {
    std::string fname = outfile_name + "/" + ref_names[i] + "_reads.txt";
    if (gzip_output) {
      outfiles[i] = std::unique_ptr<std::ostream>(new zstr::ofstream(fname + ".gz"));
    } else {
      outfiles[i] = std::unique_ptr<std::ostream>(new std::ofstream(fname));
    }
  }
  std::cout << "n_refs: " << ref_names.size() << std::endl;
  std::cout << "n_probs: " << probs.size() << std::endl;
  std::cout << "writing reads..." << std::endl;
  write_reads(reads_to_ec, probs, group_indices, outfiles);

  return 0;
}
