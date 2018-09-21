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

int main(int argc, char* argv[]) {
  std::string probs_file = std::string(GetOpt(argv, argv+argc, "-p"));
  std::string abundances_file = std::string(GetOpt(argv, argv+argc, "-a"));
  std::vector<std::string> ref_names;
  std::cout << "Reading probs" << std::endl;
  const std::unordered_map<long unsigned, std::vector<bool>> &probs = read_probs(probs_file, abundances_file, ref_names);

  std::cout << "Reading .sam file" << std::endl;
  std::string sam_file = std::string(GetOpt(argv, argv+argc, "-s"));
  const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec = read_sam(sam_file);

  std::cout << "Reading pseudoalignments" << std::endl;
  std::string ec_file = std::string(GetOpt(argv, argv+argc, "-e"));
  const std::unordered_map<std::vector<bool>, long unsigned> &ec_to_id = read_ec_ids(ec_file, reads_in_ec);

  std::cout << "Assigning reads to reference groups" << std::endl;
  std::string outfile = std::string(GetOpt(argv, argv+argc, "-o"));
  write_reads(ec_to_id, reads_in_ec, probs, ref_names, outfile);

  return 0;
}
