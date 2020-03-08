#include "mGEMS.h"

#include <algorithm>

#include "bin_reads.h"
#include "extract_bin.h"

namespace mGEMS {
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::vector<std::string> &group_names, std::istream &probs_file, std::vector<std::string> *target_groups) {
  uint32_t num_ecs = aln.size();
  uint32_t n_groups = abundances.size();
  std::vector<long double> thresholds(n_groups);
  ConstructThresholds(num_ecs, theta_frac, abundances, &thresholds);

  std::vector<std::vector<bool>> assignments(num_ecs, std::vector<bool>(n_groups, false));
  std::vector<std::vector<uint32_t>> read_bins(n_groups, std::vector<uint32_t>());
  const std::vector<bool> &groups_to_assign = AssignProbs(thresholds, probs_file, target_groups, &assignments, aln.aligned_reads, &read_bins);

  //BinReads(assignments, groups_to_assign, aln.aligned_reads, &read_bins);

  std::vector<std::vector<uint32_t>> out_bins;
  for (uint32_t i = 0; i < n_groups; ++i) {
    if (read_bins[i].size() > 0) {
      out_bins.push_back(std::move(read_bins[i]));
      std::sort(out_bins.back().begin(), out_bins.back().end());
    }
  }
  return out_bins;
}
}
