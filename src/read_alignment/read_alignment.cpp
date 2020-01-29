#include <string>
#include <iostream>

#include "cxxargs/include/cxxargs.hpp"

#include "file.hpp"
#include "read_alignment/read_files.h"
#include "read_alignment/write_files.h"

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_argument<std::string>('e', "ec-path", "pseudoalignments.ec file path.");
  args.add_short_argument<std::string>('o', "Output directory.");
  args.add_argument<std::string>('s', "sam-path", "pseudobam (.sam format) file path. (not required if --cin is specified)");
  args.set_not_required("sam-path");
  args.add_long_argument<bool>("themisto", "Themisto compatibility mode. (default: false)", false);
  args.add_long_argument<uint32_t>("n-refs", "Number of reference sequences in the pseudoalignment (only required in Themisto-mode)");
  args.set_not_required("n-refs");
  args.add_long_argument<bool>("cin", "Read the sam-file from cin. (default: false)", false);
  args.add_long_argument("gzip-output", "Compress the output files (default: false)", false);

  args.parse(argc, argv);
  if (!args.value<bool>("cin") && !args.is_initialized("sam-path")) {
    throw std::runtime_error("sam-path must be specified if not reading from cin");
  }
  if (args.value<bool>("themisto") && !args.is_initialized("n-refs")) {
    throw std::runtime_error("n-refs must be specified if running in Themisto-mode.");
  }
}

int main(int argc, char* argv[]) {
  cxxargs::Arguments args("read-alignment", "Usage:...");
  try {
    std::cerr << "Parsing arguments" << std::endl;
    parse_args(argc, argv, args);
  } catch (std::exception &e) {
    std::cerr << "Parsing arguments failed:\n"
	      << "\t" << e.what()
	      << "\nexiting\n";
    return 0;
  }
  std::cerr << "Reading alignment files" << std::endl;
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  File::In ec_file(args.value<std::string>('e'));
  if (args.value<bool>("cin")) {
    reads_in_ec(std::cin, ec_file, &reads_to_ec, args.value<bool>("themisto"), (args.value<bool>("themisto") ? args.value<uint32_t>("n-refs") : 0));
  } else {
    std::cerr << "reading from: "  << args.value<std::string>('s') << std::endl;
    File::In sam_file(args.value<std::string>('s'));
    reads_in_ec(sam_file.stream(), ec_file, &reads_to_ec, args.value<bool>("themisto"), (args.value<bool>("themisto") ? args.value<uint32_t>("n-refs") : 0));
  }

  std::cerr << "Writing read assignments to equivalence classes" << std::endl;
  File::Out outfile;
  if (args.value<bool>("gzip-output")) {
    outfile.open_compressed(args.value<std::string>('o') + "/" + "ec_to_read.csv" + ".gz");
  } else {
    outfile.open(args.value<std::string>('o') + "/" + "ec_to_read.csv");
  }
  write_ecs(reads_to_ec, outfile);

  return 0;
}
