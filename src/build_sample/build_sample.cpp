#include <string>
#include <algorithm>

#include "arguments.h"
#include "build_sample/assign_reads.h"
#include "zstr/zstr.hpp"

int main(int argc, char* argv[]) {
  std::string assignment_path = std::string(GetOpt(argv, argv+argc, "-a"));
  std::string outfile = std::string(GetOpt(argv, argv+argc, "-o"));
  std::string strand1 = std::string(GetOpt(argv, argv+argc, "-1"));
  std::string strand2 = std::string(GetOpt(argv, argv+argc, "-2"));
  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");

  zstr::ifstream assignment_file(assignment_path);
  const std::map<std::string, std::set<short unsigned>> &assignments = read_assignments(assignment_file, 0);

  std::unique_ptr<std::istream> infiles[2];
  infiles[0] = std::unique_ptr<std::istream>(new zstr::ifstream(strand1));
  infiles[1] = std::unique_ptr<std::istream>(new zstr::ifstream(strand2));
  std::unique_ptr<std::ostream> outfiles[2];
  if (gzip_output) {
    outfiles[0] = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile + "_1.fastq.gz"));
    outfiles[1] = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile + "_2.fastq.gz"));
  } else {
    outfiles[0] = std::unique_ptr<std::ostream>(new std::ofstream(outfile + "_1.fastq"));
    outfiles[1] = std::unique_ptr<std::ostream>(new std::ofstream(outfile + "_2.fastq"));
  }

  assign_reads(assignments, outfiles, infiles);
  
  return 0;
}
