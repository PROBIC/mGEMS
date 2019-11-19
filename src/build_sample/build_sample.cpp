#include <string>
#include <algorithm>
#include <exception>
#include <iostream>

#include "cxxargs/include/cxxargs.hpp"
#include "zstr/zstr.hpp"

#include "arguments.h"
#include "build_sample/assign_reads.h"

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('a', "Read assignments to equivalence classes.");
  args.add_short_argument<std::string>('o', "Output directory.");
  args.add_short_argument<std::string>('1', "Strand 1 to demix reads from.");
  args.add_short_argument<std::string>('2', "Strand 2 to demix reads from.");
  args.add_long_argument<bool>("gzip-output", "Compress output files (default: false)");
  args.parse(argc, argv);
}

int main(int argc, char* argv[]) {
  cxxargs::Arguments args("build-sample", "Usage:...");
  try {
    std::cout << "Parsing arguments" << std::endl;
    parse_args(argc, argv, args);
  } catch (std::exception &e) {
    std::cout << "Parsing arguments failed:\n"
	      << '\t' << e.what()
	      << "\nexiting\n";
    return 0;
  }

  zstr::ifstream assignment_file(args.value<std::string>('a'));
  const std::set<std::string> &assignments = read_assignments(assignment_file);

  std::unique_ptr<std::istream> infiles[2];
  infiles[0] = std::unique_ptr<std::istream>(new zstr::ifstream(args.value<std::string>('1')));
  infiles[1] = std::unique_ptr<std::istream>(new zstr::ifstream(args.value<std::string>('2')));
  std::unique_ptr<std::ostream> outfiles[2];
  std::string outfile = args.value<std::string>('o');
  if (args.value<bool>("gzip-output")) {
    outfiles[0] = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile + "_1.fastq.gz"));
    outfiles[1] = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile + "_2.fastq.gz"));
  } else {
    outfiles[0] = std::unique_ptr<std::ostream>(new std::ofstream(outfile + "_1.fastq"));
    outfiles[1] = std::unique_ptr<std::ostream>(new std::ofstream(outfile + "_2.fastq"));
  }

  assign_reads(assignments, outfiles, infiles);
  
  return 0;
}
