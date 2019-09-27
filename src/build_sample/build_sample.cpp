#include <string>
#include <algorithm>

#include "build_sample/assign_reads.h"
#include "zstr/zstr.hpp"

char* GetOpt(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

int main(int argc, char* argv[]) {
  std::string assignment_file = std::string(GetOpt(argv, argv+argc, "-a"));
  std::string outfile = std::string(GetOpt(argv, argv+argc, "-o"));
  std::string strand1 = std::string(GetOpt(argv, argv+argc, "-1"));
  std::string strand2 = std::string(GetOpt(argv, argv+argc, "-2"));
  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");

  const std::map<std::string, std::set<short unsigned>> &assignments = read_assignments(assignment_file, 0);
  
  assign_reads(outfile, strand1, strand2, gzip_output, assignments);
  
  return 0;
}
