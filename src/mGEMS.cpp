#include <exception>
#include <iostream>
#include <cstring>

#include "telescope.hpp"
#include "cxxargs.hpp"
#include "cxxio.hpp"

#include "bin_reads.h"
#include "extract_bin.h"

void ParseExtract(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::vector<std::string>>('r', "Input reads (comma separated list)");
  args.add_short_argument<std::string>('o', "Output directory.");
  args.add_long_argument<std::vector<std::string>>("bins", "Comma-separated list of bins to extract from the paired-end reads.");
  args.add_long_argument<bool>("compress", "Compress extracted reads with zlib (.gz extension, default: true)", true);
  args.add_long_argument<bool>("write-unassigned", "Extract reads that pseudoaligned to a reference sequence but were not assigned to any group.", false);
  args.add_long_argument<bool>("write-assignment-table", "Write the read-group assignment table.", false);
  args.set_not_required("bins");
  args.set_not_required('o');

  args.parse(argc, argv);
}

void ParseBin(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_long_argument<std::vector<std::string>>("themisto-alns", "Comma-separated list of paired-end alignments from Themisto.");
  args.add_short_argument<std::string>('o', "Output directory (must exist before running).");
  args.add_short_argument<std::string>('a', "Relative abundances estimates from mSWEEP.");
  args.add_long_argument<std::string>("probs", "Posterior probabilities from mSWEEP.");
  args.add_long_argument<std::string>("merge-mode", "How to merge paired-end alignments from Themisto (default: intersection).", "intersection");
  args.add_long_argument<std::vector<std::string>>("groups", "Which reference groups to bin reads to (default: all).");
  args.add_long_argument<long double>("min-abundance", "Bin only the groups that have a relative abundance higher than this value (optional).");
  args.add_short_argument<long double>('q', "Tuning parameter for the binning thresholds (default: 1.0).", (long double)1);
  args.add_long_argument<std::string>("index", "Themisto pseudoalignment index directory.");
  args.add_long_argument<bool>("write-unassigned", "Extract reads that pseudoaligned to a reference sequence but were not assigned to any group.", false);
  args.add_long_argument<bool>("write-assignment-table", "Write the read-group assignment table.", false);
  args.set_not_required("groups");
  args.set_not_required("min-abundance");

  args.parse(argc, argv);
}

void Extract(const std::vector<std::vector<uint32_t>> &bins, const std::vector<uint32_t> &unassigned_bin, const std::vector<std::string> &target_groups, const cxxargs::Arguments &args) {
  uint32_t n_out_groups = bins.size();
  uint8_t n_strands = args.value<std::vector<std::string>>('r').size();
  std::vector<cxxio::In> in_strands(n_strands);
  std::vector<cxxio::Out> out_strands(n_strands);
  for (uint32_t i = 0; i < n_out_groups; ++i) {
    for (uint8_t j = 0; j < n_strands; ++j) {
      in_strands[j].open(args.value<std::vector<std::string>>('r')[j]);
      std::string out_name = args.value<std::string>('o') + "/" + target_groups[i] + "_" + std::to_string(j + 1) + ".fastq";
      if (args.value<bool>("compress")) {
	out_name += ".gz";
	out_strands[j].open_compressed(out_name);
      } else {
	out_strands[j].open(out_name);
      }
    }
    if (bins[i].size() > 0) {
      mGEMS::ExtractBin(bins[i], in_strands, &out_strands);
    } else {
      std::cerr << "WARNING: No reads assigned to bin " << target_groups[i] << '.' << std::endl;
    }
  }
  if (args.value<bool>("write-unassigned")) {
    // Extract unassigned reads
    for (uint8_t j = 0; j < n_strands; ++j) {
      in_strands[j].open(args.value<std::vector<std::string>>('r')[j]);
      std::string out_name = args.value<std::string>('o') + '/' + "unassigned_reads" + '_' + std::to_string(j + 1) + ".fastq";
      if (args.value<bool>("compress")) {
	out_name += ".gz";
	out_strands[j].open_compressed(out_name);
      } else {
	out_strands[j].open(out_name);
      }
      if (unassigned_bin.size() > 0) {
	mGEMS::ExtractBin(unassigned_bin, in_strands, &out_strands);
      }
    }
  }
}

void ReadAndExtract(cxxargs::Arguments &args) {
  cxxio::directory_exists(args.value<std::string>('o'));
  uint32_t n_bins = args.value<std::vector<std::string>>("bins").size();
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::string> target_groups(n_bins);
  for (uint32_t i = 0; i < n_bins; ++i) {
    cxxio::In istream(args.value<std::vector<std::string>>("bins")[i]);
    bins.emplace_back(mGEMS::ReadBin(istream.stream()));
    std::string out_name = args.value<std::vector<std::string>>("bins")[i];
    if (out_name.find(".") != std::string::npos) {
      out_name.erase(out_name.find("."), out_name.size());
    }
    if (out_name.rfind("/") != std::string::npos) {
      out_name.erase(0, out_name.rfind("/") + 1);
    }
    target_groups[i] = out_name;
  }
  Extract(bins, std::vector<uint32_t>(), target_groups, args);
}

void FilterTargetGroups(const std::vector<std::string> &group_names, const std::vector<long double> &abundances, const long double min_abundance, std::vector<std::string> *target_groups) {
  uint32_t n_groups = group_names.size();
  for (uint32_t i = 0; i < n_groups; ++i) {
    if (abundances[i] < min_abundance && std::find(target_groups->begin(), target_groups->end(), group_names[i]) != target_groups->end()) {
      target_groups->erase(std::find(target_groups->begin(), target_groups->end(), group_names[i]));
    }
  }
}

void Bin(const cxxargs::Arguments &args, bool extract_bins) {
  cxxio::directory_exists(args.value<std::string>('o').c_str());
  cxxio::directory_exists(args.value<std::string>("index").c_str());

  std::vector<std::string> groups;
  std::vector<long double> abundances;
  cxxio::In msweep_abundances(args.value<std::string>('a'));
  mGEMS::ReadAbundances(msweep_abundances.stream(), &abundances, &groups);
  
  std::vector<std::istream*> themisto_alns;
  for (uint32_t i = 0; i < args.value<std::vector<std::string>>("themisto-alns").size(); ++i) {
    cxxio::In(args.value<std::vector<std::string>>("themisto-alns")[i]);
    themisto_alns.emplace_back(new bxz::ifstream(args.value<std::vector<std::string>>("themisto-alns")[i]));
  }

  cxxio::In themisto_index(args.value<std::string>("index") + "/coloring-names.txt");
  uint32_t n_refs = themisto_index.count_lines<uint32_t>();
  themisto_index.close();
  ThemistoAlignment aln;
  ReadThemisto(get_mode(args.value<std::string>("merge-mode")), n_refs, themisto_alns, &aln);

  cxxio::In probs_file(args.value<std::string>("probs"));
  std::vector<std::string> target_groups;
  if (args.is_initialized("groups")) {
    target_groups = args.value<std::vector<std::string>>("groups");
  } else {
    target_groups = groups;
  }
  if (args.is_initialized("min-abundance")) {
    FilterTargetGroups(groups, abundances, args.value<long double>("min-abundance"), &target_groups);
  }

  uint32_t n_groups = abundances.size();
  std::vector<uint32_t> unassigned_bin;
  std::vector<std::vector<bool>> assignments_mat(aln.size(), std::vector<bool>(n_groups, false));
  const std::vector<std::vector<uint32_t>> &bins = mGEMS::Bin(aln, args.value<long double>('q'), abundances, probs_file.stream(), &target_groups, &unassigned_bin, &assignments_mat);

  std::cerr << args.value<bool>("write-assignment-table") << std::endl;
  if (args.value<bool>("write-assignment-table")) {
    cxxio::Out of(args.value<std::string>('o') + '/' + "reads_to_groups.tsv");
    of.stream() << "#read" << '\t';
    for (uint32_t i = 0; i < n_groups; ++i) {
      of.stream() << groups[i];
      if (i < n_groups - 1) {
	of.stream() << '\t';
      }
    }
    of.stream() << '\n';
    mGEMS::WriteAssignments(assignments_mat, aln, of.stream());
    of.close();
  }

  if (!extract_bins) {
    uint32_t n_bins = bins.size();
    for (uint32_t i = 0; i < n_bins; ++i) {
      cxxio::Out of(args.value<std::string>('o') + '/' + target_groups[i] + ".bin");
      mGEMS::WriteBin(bins[i], of.stream());
    }
    if (args.value<bool>("write-unassigned")) {
      cxxio::Out of(args.value<std::string>('o') + '/' + "unassigned_reads.bin");
      mGEMS::WriteBin(unassigned_bin, of.stream());
    }
  } else {
    Extract(bins, unassigned_bin, target_groups, args);
  }
}

int main(int argc, char* argv[]) {
  cxxargs::Arguments args("mGEMS", "Usage: mGEMS -r <input-reads_1>,<input-reads_2> --themisto-alns <input-reads_1 pseudoalignments>,<input-reads_2 pseudoalignments> -o <output directory> --probs <posterior probabilities> -a <abundance estimates> --index <Themisto index> --groups <group names to extract (optional)>");
  try {
    if (argc < 2) {
      std::cerr << args.help() << std::endl;
      return 0;
    } else if (strcmp(argv[1], "bin") == 0) {
      ParseBin(argc, argv, args);
      Bin(args, false);
    } else if (strcmp(argv[1], "extract") == 0) {
      ParseExtract(argc, argv, args);
      if (!args.is_initialized("bins")) {
	throw std::runtime_error("Bins must be specified.");
      }
      if (!args.is_initialized('o')) {
	throw std::runtime_error("Output directory not given.");
      }
      ReadAndExtract(args);
    } else {
      ParseBin(argc, argv, args);
      ParseExtract(argc, argv, args);
      Bin(args, true);
    }
  } catch (std::exception &e) {
    std::cerr << "Error:"
	      << "\t" << e.what()
	      << "\nexiting\n";
    return 1;
  }

  return 0;
}
