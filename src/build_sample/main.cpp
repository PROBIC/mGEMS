#include <string>
#include <algorithm>

#include "build_sample/assign_reads.h"

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

  assign_reads(assignment_file, outfile, strand1, strand2);
  
  return 0;
}
