#include "bin_reads.h"

#include <sstream>
#include <iterator>
#include <algorithm>
#include <utility>
#include <cmath>

namespace mGEMS {
uint32_t ReadAbundances(std::istream &stream, std::vector<long double> *abundances, std::vector<std::string> *groups) {
  std::string line;
  uint32_t n_groups = 0;
  while(std::getline(stream, line)) {
    if (line.at(0) != '#') {
      std::string part;
      std::stringstream parts(line);
      std::getline(parts, part, '\t');
      groups->emplace_back(part);
      std::getline(parts, part, '\t');
      abundances->emplace_back(std::stold(part));
      ++n_groups;
    }
  }
  return n_groups;
}

void ConstructThresholds(const uint32_t num_ecs, const long double theta_frac, const std::vector<long double> &abundances, std::vector<long double> *thresholds) {
  long double log_thresh = std::log1pl(-(long double)abundances.size()/(long double)num_ecs);
  log_thresh += std::log(theta_frac);
  uint32_t n_refs = abundances.size();
  for (uint32_t i = 0; i < n_refs; ++i) {
    (*thresholds)[i] = std::exp(std::log(abundances[i]) + log_thresh);
  }
}

void MaskProbs(const std::string &header_line, std::vector<std::string> *target_groups, std::vector<bool> *mask) {
  std::string part;
  std::stringstream partition(header_line);
  std::getline(partition, part, ',');
  uint32_t ref_id = 0;
  std::vector<std::string> ordered_targets(target_groups->size());
  uint32_t target_id = 0;
  while (std::getline(partition, part, ',')) {
    if (std::find(target_groups->begin(), target_groups->end(), part) != target_groups->end()) {
      (*mask)[ref_id] = true;
      ordered_targets[target_id] = part;
      ++target_id;
    }
    ++ref_id;
  }
  (*target_groups) = ordered_targets;
}

void AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, const std::vector<bool> &mask, std::vector<std::vector<bool>> *assignments, const std::vector<std::vector<uint32_t>> &aligned_reads, std::vector<std::vector<uint32_t>> *bins) {
  std::vector<uint32_t> n_reads(thresholds.size(), 0);
  uint32_t ec_id = 0;

  std::string line;
  while (std::getline(probs_file, line)) {
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ','); // first element is the ec id
    uint32_t ref_id = 0;
    while(std::getline(partition, part, ',')) {
      long double abundance = std::stold(part);
      (*assignments)[ec_id][ref_id] = (abundance >= thresholds[ref_id]) && mask[ref_id];
      if ((*assignments)[ec_id][ref_id]) {
	for (uint32_t i = 0; i < aligned_reads[ec_id].size(); ++i) {
	  (*bins)[ref_id].push_back(aligned_reads[ec_id][i] + 1);
	}
      }
      ++ref_id;
    }
    ++ec_id;
  }
}

void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of) {
  for (uint32_t i = 0; i < binned_reads.size(); ++i) {
    of << binned_reads[i] << '\n';
  }
  of.flush();
}

std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::istream &probs_file, std::vector<std::string> *target_groups) {
  uint32_t num_ecs = aln.size();
  uint32_t n_groups = abundances.size();
  std::vector<long double> thresholds(n_groups);
  ConstructThresholds(num_ecs, theta_frac, abundances, &thresholds);

  std::vector<std::vector<bool>> assignments(num_ecs, std::vector<bool>(n_groups, false));
  std::vector<std::vector<uint32_t>> read_bins(n_groups, std::vector<uint32_t>());

  std::string probs_header;
  std::getline(probs_file, probs_header); // 1st line is header
  std::vector<bool> mask(n_groups, false);
  MaskProbs(probs_header, target_groups, &mask);

  AssignProbs(thresholds, probs_file, mask, &assignments, aln.aligned_reads, &read_bins);

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
