#include <string>
#include <iostream>

#include "arguments.h"
#include "zstr/zstr.hpp"
#include "read_alignment/read_files.hpp"
#include "read_alignment/write_files.hpp"

int main(int argc, char* argv[]) {
  std::cout << "Reading alignment files" << std::endl;
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  std::string ec_path = std::string(GetOpt(argv, argv+argc, "-e"));
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  if (read_from_cin) {
    reads_in_ec(std::cin, ec_path, &reads_to_ec);
  } else {
    std::string sam_path = std::string(GetOpt(argv, argv+argc, "-s"));
    zstr::ifstream sam_file(sam_path);
    reads_in_ec(sam_file, ec_path, &reads_to_ec);
  }

  std::cout << "Writing read assignments to equivalence classes" << std::endl;
  std::string outfile_name = std::string(GetOpt(argv, argv+argc, "-o"));
  std::unique_ptr<std::ostream> outfile;
  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");
  if (gzip_output) {
    outfile = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile_name + "/" + "ec_to_read.csv" + ".gz"));
  } else {
    outfile = std::unique_ptr<std::ostream>(new std::ofstream(outfile_name + "/" + "ec_to_read.csv"));
  }
  write_ecs(reads_to_ec, outfile);

  return 0;
}
