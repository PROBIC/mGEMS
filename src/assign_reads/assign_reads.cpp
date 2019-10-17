#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>

#include "arguments.h"
#include "assign_reads/read_files.h"
#include "assign_reads/write_files.h"
#include "assign_reads/construct_thresholds.h"
#include "zstr/zstr.hpp"

int main(int argc, char* argv[]) {
  std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> reads_in_ec;
  std::cout << "Reading read assignments to equivalence classes" << std::endl;
  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::unique_ptr<std::istream> assignments_file = OpenInstream(argv, argv+argc, "-f");
    read_ec_assignments(*assignments_file, &reads_in_ec);
  }

  std::cout << "Constructing assignment thresholds from abundances" << std::endl;
  std::vector<std::pair<std::string, long double>> thresholds;
  if (CmdOptionPresent(argv, argv+argc, "-a")) {
    std::unique_ptr<std::istream> abundances_file = OpenInstream(argv, argv+argc, "-a");
    long double mult = (CmdOptionPresent(argv, argv+argc, "-q") ? std::stold(std::string(GetOpt(argv, argv+argc, "-q"))) : (long double)1.0);
    construct_thresholds(reads_in_ec.size(), *abundances_file, &thresholds, mult);
  }
  uint16_t n_refs = thresholds.size();

  std::cout << "Reading probs" << std::endl;
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  if (read_from_cin) {
    assign_reads(thresholds, std::cin, &reads_in_ec);
  } else {
    std::unique_ptr<std::istream> probs_file = OpenInstream(argv, argv+argc, "-p");
    assign_reads(thresholds, *probs_file, &reads_in_ec);
  }

  std::cout << "Assigning reads to reference groups" << std::endl;
  std::vector<short unsigned> group_indices;
  std::unique_ptr<std::istream> groups_file = OpenInstream(argv, argv+argc, "--groups");
  read_groups_filter(thresholds, *groups_file, &group_indices);
  
  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");
  std::string outfile_name = std::string(GetOpt(argv, argv+argc, "-o"));
  std::vector<std::unique_ptr<std::ostream>> outfiles(n_refs);
  for (auto i : group_indices) {
    std::string fname = outfile_name + "/" + thresholds.at(i).first + "_reads.txt";
    if (gzip_output) {
      outfiles[i] = std::unique_ptr<std::ostream>(new zstr::ofstream(fname + ".gz"));
    } else {
      outfiles[i] = std::unique_ptr<std::ostream>(new std::ofstream(fname));
    }
  }
  std::cout << "writing reads..." << std::endl;
  write_reads(reads_in_ec, group_indices, outfiles);

  return 0;
}
