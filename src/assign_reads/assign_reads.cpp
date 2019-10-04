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

void multiply_abundances(std::vector<std::pair<std::string, long double>> &abundances, long double log_thresh) {
  for (size_t i = 0; i < abundances.size(); ++i) {
    abundances.at(i).second = std::exp(std::log(abundances.at(i).second) + log_thresh);
  }
}

int main(int argc, char* argv[]) {
  std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> reads_to_ec;
  std::cout << "Reading read assignments to equivalence classes" << std::endl;
  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::unique_ptr<std::istream> assignments_file = OpenInstream(argv, argv+argc, "-f");
    read_assignments(*assignments_file, &reads_to_ec);
  }

  std::cout << "Reading abundances" << std::endl;
  std::unique_ptr<std::istream> abundances_file = OpenInstream(argv, argv+argc, "-a");
  std::vector<std::pair<std::string, long double>> abundances = read_abundances(*abundances_file);
  uint16_t n_refs = abundances.size();
  double log_thresh = std::log1pl(-(long double)n_refs/(long double)reads_to_ec.size());
  if (CmdOptionPresent(argv, argv+argc, "-q")) {
    log_thresh += std::stold(std::string(GetOpt(argv, argv+argc, "-q")));
  }
  std::cout << std::exp(log_thresh) << std::endl;
  multiply_abundances(abundances, log_thresh);
    
  std::cout << "Reading probs" << std::endl;
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  if (read_from_cin) {
    read_probs(abundances, std::cin, &reads_to_ec);
  } else {
    std::unique_ptr<std::istream> probs_file = OpenInstream(argv, argv+argc, "-p");
    read_probs(abundances, *probs_file, &reads_to_ec);
  }

  bool all_groups = !CmdOptionPresent(argv, argv+argc, "--groups");
  std::cout << "Assigning reads to reference groups" << std::endl;
  std::vector<short unsigned> group_indices;
  if (all_groups) {
    group_indices.resize(n_refs);
    for (size_t i = 0; i < n_refs; ++i) {
      group_indices[i] = i;
    }
  } else {
    std::unique_ptr<std::istream> groups_file = OpenInstream(argv, argv+argc, "--groups");
    read_groups(abundances, *groups_file, &group_indices);
  }
  
  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");
  std::string outfile_name = std::string(GetOpt(argv, argv+argc, "-o"));
  std::vector<std::unique_ptr<std::ostream>> outfiles(n_refs);
  for (auto i : group_indices) {
    std::string fname = outfile_name + "/" + abundances.at(i).first + "_reads.txt";
    if (gzip_output) {
      outfiles[i] = std::unique_ptr<std::ostream>(new zstr::ofstream(fname + ".gz"));
    } else {
      outfiles[i] = std::unique_ptr<std::ostream>(new std::ofstream(fname));
    }
  }
  std::cout << "writing reads..." << std::endl;
  write_reads(reads_to_ec, group_indices, outfiles);

  return 0;
}
