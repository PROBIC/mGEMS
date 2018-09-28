#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>

#include "read_files.h"
#include "write_files.h"

char* GetOpt(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

int main(int argc, char* argv[]) {
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  std::cout << "Reading read assignments to equivalence classes" << std::endl;
  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::string assignments_file = std::string(GetOpt(argv, argv+argc, "-f"));
    reads_to_ec = read_assignments(assignments_file);
  } else {
    std::string sam_file = std::string(GetOpt(argv, argv+argc, "-s"));
    std::string ec_file = std::string(GetOpt(argv, argv+argc, "-e"));
    reads_to_ec = reads_in_ec(sam_file, ec_file);
  }

  std::string outfile = std::string(GetOpt(argv, argv+argc, "-o"));
  if (CmdOptionPresent(argv, argv+argc, "--write-ecs")) {
    std::cout << "Writing read assignments to equivalence classes" << std::endl;
    write_ecs(reads_to_ec, outfile);
  } else {
    double theta_frac = 1.0;
    if (CmdOptionPresent(argv, argv+argc, "-q")) {
      theta_frac = std::stod(std::string(GetOpt(argv, argv+argc, "-q")));
    }
    std::string probs_file = std::string(GetOpt(argv, argv+argc, "-p"));
    std::string abundances_file = std::string(GetOpt(argv, argv+argc, "-a"));
    std::vector<std::string> ref_names;
    std::cout << "Reading probs" << std::endl;
    const std::unordered_map<long unsigned, std::vector<bool>> &probs = read_probs(probs_file, abundances_file, ref_names, theta_frac);

    bool all_groups = !CmdOptionPresent(argv, argv+argc, "--groups");

    std::cout << "Assigning reads to reference groups" << std::endl;
    if (all_groups) {
      std::vector<short unsigned> group_indices(ref_names.size());
      for (size_t i = 0; i < ref_names.size(); ++i) {
	group_indices[i] = i;
      }
      write_reads(reads_to_ec, probs, ref_names, group_indices, outfile);
    } else {
      std::string groups_file = std::string(GetOpt(argv, argv+argc, "--groups"));
      std::vector<short unsigned> group_indices = read_groups(groups_file, ref_names);
      write_reads(reads_to_ec, probs, ref_names, group_indices, outfile);
    }
  }

  return 0;
}
