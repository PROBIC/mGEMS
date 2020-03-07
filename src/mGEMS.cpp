#include "mGEMS.h"

#include "assign_reads.h"
#include "extract_bin.h"

namespace mGEMS {
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::vector<std::string> &group_names, std::istream &probs_file, std::vector<std::string> *target_groups) {
  uint32_t num_ecs = aln.size();
  uint32_t n_groups = abundances.size();
  std::vector<long double> thresholds(n_groups);
  ConstructThresholds(num_ecs, theta_frac, abundances, &thresholds);

  std::vector<std::vector<bool>> assignments(num_ecs, std::vector<bool>(n_groups, false));
  AssignProbs(thresholds, probs_file, &assignments);
  
  uint32_t n_out_groups = target_groups->size();
  std::vector<bool> groups_to_assign(n_groups, false);
  std::vector<std::string> order_in_output(n_out_groups);
  for (uint32_t i = 0; i < n_out_groups; ++i) {
    std::vector<std::string>::iterator it = std::find(group_names.begin(), group_names.end(), (*target_groups)[i]);
    int index = std::distance(group_names.begin(), it);
    groups_to_assign[index] = true;
    order_in_output[i] = group_names[index];
  }
  (*target_groups) = order_in_output;

  std::vector<std::vector<uint32_t>> read_bins(n_out_groups);
  BinReads(assignments, groups_to_assign, aln.aligned_reads, &read_bins);

  for (uint32_t i = 0; i < n_out_groups; ++i) {
    std::sort(read_bins[i].begin(), read_bins[i].end());
  }
  return read_bins;
}

void Extract(const std::vector<std::string> &target_groups, const std::string &outdir, const std::string &strand_1, const std::string &strand_2, std::vector<std::vector<uint32_t>> &bins) {
  uint32_t n_out_groups = target_groups.size();
  for (uint32_t i = 0; i < n_out_groups; ++i) {
    File::In istrand_1(strand_1);
    File::In istrand_2(strand_2);
    File::Out ostrand_1(outdir + target_groups[i] + "_1.fastq");
    File::Out ostrand_2(outdir + target_groups[i] + "_2.fastq");
    std::sort(bins[i].begin(), bins[i].end());
    ExtractBin(bins[i], &ostrand_1.stream(), &ostrand_2.stream(), &istrand_1.stream(), &istrand_2.stream());
  }
}
}
