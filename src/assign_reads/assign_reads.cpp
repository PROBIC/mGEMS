#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>

#include "assign_reads/read_files.h"
#include "assign_reads/write_files.h"
#include "zstr/zstr.hpp"

char* GetOpt(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

int main(int argc, char* argv[]) {
  std::unordered_map<long unsigned, std::vector<std::string>> reads_to_ec;
  bool read_from_cin = CmdOptionPresent(argv, argv+argc, "--cin");
  std::cout << "Reading read assignments to equivalence classes" << std::endl;
  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::string assignments_path = std::string(GetOpt(argv, argv+argc, "-f"));
    std::ifstream assignments_file(assignments_path);
    read_assignments(assignments_file, &reads_to_ec);
    assignments_file.close();
  } else {
    std::string sam_path = std::string(GetOpt(argv, argv+argc, "-s"));
    std::string ec_path = std::string(GetOpt(argv, argv+argc, "-e"));
    zstr::ifstream ec_file(ec_path);
    if (read_from_cin) {
      reads_in_ec(std::cin, ec_file, &reads_to_ec);
    } else {
      zstr::ifstream sam_file(sam_path);
      reads_in_ec(sam_file, ec_file, &reads_to_ec);
    }
  }

  bool gzip_output = CmdOptionPresent(argv, argv+argc, "--gzip-output");

  std::string outfile_name = std::string(GetOpt(argv, argv+argc, "-o"));
  if (CmdOptionPresent(argv, argv+argc, "--write-ecs")) {
    std::cout << "Writing read assignments to equivalence classes" << std::endl;
    std::unique_ptr<std::ostream> outfile;
    if (gzip_output) {
      outfile = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile_name + "/" + "ec_to_read.csv" + ".gz"));
    } else {
      outfile = std::unique_ptr<std::ostream>(new std::ofstream(outfile_name + "/" + "ec_to_read.csv"));
    }
    write_ecs(reads_to_ec, outfile);
  } else {
    double theta_frac = 1.0;
    if (CmdOptionPresent(argv, argv+argc, "-q")) {
      theta_frac = std::stod(std::string(GetOpt(argv, argv+argc, "-q")));
    }
    std::string probs_file = std::string(GetOpt(argv, argv+argc, "-p"));
    std::string abundances_file = std::string(GetOpt(argv, argv+argc, "-a"));
    std::vector<std::string> ref_names;
    std::cout << "Reading probs" << std::endl;
    const std::unordered_map<long unsigned, std::vector<bool>> &probs = read_probs(probs_file, abundances_file, ref_names, theta_frac);

    bool all_groups = !CmdOptionPresent(argv, argv+argc, "--groups");

    std::cout << "Assigning reads to reference groups" << std::endl;
    std::vector<short unsigned> group_indices;
    if (all_groups) {
      group_indices.resize(ref_names.size());
      for (size_t i = 0; i < ref_names.size(); ++i) {
	group_indices[i] = i;
      }
    } else {
      std::string groups_path = std::string(GetOpt(argv, argv+argc, "--groups"));
      std::ifstream groups_file(groups_path);
      read_groups(ref_names, groups_file, &group_indices);
    }
    std::vector<std::unique_ptr<std::ostream>> outfiles(ref_names.size());
    for (auto i : group_indices) {
      std::string fname = outfile_name + "/" + ref_names[i] + "_reads.txt";
      if (gzip_output) {
	outfiles[i] = std::unique_ptr<std::ostream>(new zstr::ofstream(fname + ".gz"));
      } else {
	outfiles[i] = std::unique_ptr<std::ostream>(new std::ofstream(fname));
      }
    }
    std::cout << "n_refs: " << ref_names.size() << std::endl;
    std::cout << "n_probs: " << probs.size() << std::endl;
    std::cout << "writing reads..." << std::endl;
    write_reads(reads_to_ec, probs, group_indices, outfiles);
  }

  return 0;
}
