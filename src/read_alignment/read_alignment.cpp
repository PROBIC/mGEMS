#include <string>
#include <iostream>

#include "arguments.h"
#include "zstr/zstr.hpp"
#include "read_alignment/read_files.h"
#include "read_alignment/write_files.h"

int main(int argc, char* argv[]) {
  std::cout << "Reading alignment files" << std::endl;
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  std::string ec_path = std::string(GetOpt(argv, argv+argc, "-e"));
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  bool themisto = CmdOptionPresent(argv, argv+argc, "--themisto");
  if (read_from_cin) {
    reads_in_ec(std::cin, ec_path, &reads_to_ec, themisto, (themisto ? std::stoi(GetOpt(argv, argv+argc, "--n-refs")) : 0));
  } else {
    std::string sam_path = std::string(GetOpt(argv, argv+argc, "-s"));
    std::cout << "reading from: "  << sam_path << std::endl;
    zstr::ifstream sam_file(sam_path);
    std::cout << "which is: " << sam_file.good() << std::endl;
    reads_in_ec(sam_file, ec_path, &reads_to_ec, themisto, (themisto ? std::stoi(GetOpt(argv, argv+argc, "--n-refs")) : 0));
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
