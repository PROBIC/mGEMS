#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <exception>

#include "zstr/src/zstr.hpp"
#include "cxxargs/include/cxxargs.hpp"

#include "file.hpp"
#include "assign_reads/read_files.h"
#include "assign_reads/write_files.h"
#include "assign_reads/construct_thresholds.h"

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('f', "Read-EC assignments file path.");
  args.add_short_argument<std::string>('o', "Output directory.");
  args.add_short_argument<std::string>('a', "mSWEEP relative abundance estimates path.");
  args.add_short_argument<long double>('q', "Multiplier term in assignment thresholds. (default: 1.0)", (long double)1.0);
  args.add_long_argument<std::string>("groups", "Reference groups to assign reads to.");
  args.add_short_argument<std::string>('p', "Read-EC probability matrix (not required if --cin is specified)");
  args.set_not_required('p');
  args.add_long_argument<bool>("cin", "Read the probs file from cin. (default: false)", false);
  args.add_long_argument("gzip-output", "Compress the output files (default: false)", false);

  args.parse(argc, argv);
  if (!args.value<bool>("cin") && !args.is_initialized('p')) {
    throw std::runtime_error("The Read-EC probability matrix must be specified if not read from cin");
  }
}

int main(int argc, char* argv[]) {
  cxxargs::Arguments args("assign-reads", "Usage:...");
  try {
    std::cerr << "Parsing arguments" << std::endl;
    parse_args(argc, argv, args);
  } catch (std::exception &e) {
    std::cerr << "Parsing arguments failed:\n"
	      << '\t' << e.what()
	      << "\nexiting\n";
    return 0;    
  }

  std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> reads_in_ec;
  std::cerr << "Reading read assignments to equivalence classes" << std::endl;
  File::In assignments_file(args.value<std::string>('f'));
  read_ec_assignments(assignments_file.stream(), &reads_in_ec);

  std::cerr << "Constructing assignment thresholds from abundances" << std::endl;
  std::vector<std::pair<std::string, long double>> thresholds;
  File::In abundances_file(args.value<std::string>('a'));
  construct_thresholds(reads_in_ec.size(), abundances_file.stream(), &thresholds, args.value<long double>('q'));
  uint16_t n_refs = thresholds.size();

  std::cerr << "Reading probs" << std::endl;
  if (args.value<bool>("cin")) {
    assign_reads(thresholds, std::cin, &reads_in_ec);
  } else {
    File::In probs_file(args.value<std::string>('p'));
    assign_reads(thresholds, probs_file.stream(), &reads_in_ec);
  }

  std::cerr << "Assigning reads to reference groups" << std::endl;
  std::vector<short unsigned> group_indices;
  File::In groups_file(args.value<std::string>("groups"));
  read_groups_filter(thresholds, groups_file.stream(), &group_indices);
  
  std::string outfile_name = args.value<std::string>('o');
  std::vector<File::Out> outfiles(n_refs);
  for (auto i : group_indices) {
    std::string fname = outfile_name + "/" + thresholds.at(i).first + "_reads.txt";
    if (args.value<bool>("gzip-output")) {
      outfiles[i].open_compressed(fname + ".gz");
    } else {
      outfiles[i].open(fname);
    }
  }
  std::cerr << "writing reads..." << std::endl;
  write_reads(reads_in_ec, group_indices, outfiles);

  return 0;
}
