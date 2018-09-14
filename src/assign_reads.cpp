#include <string>
#include <algorithm>

#include "read_files.h"

char* GetOpt(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

int main(int argc, char* argv[]) {
  std::string sam_file = std::string(GetOpt(argv, argv+argc, "-s"));
  std::string ec_file = std::string(GetOpt(argv, argv+argc, "-p"));
  std::string outfile = std::string(GetOpt(argv, argv+argc, "-o"));
  
  read_sam(sam_file, ec_file, outfile);
  return 0;
}
