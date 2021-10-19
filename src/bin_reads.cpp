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

std::vector<bool> AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, std::vector<std::string> *target_groups, std::vector<std::vector<bool>> *assignments, const std::vector<std::vector<uint32_t>> &assigned_reads, std::vector<std::vector<uint32_t>> *bins) {
  std::string line;
  std::getline(probs_file, line); // 1st line is header
  std::vector<bool> mask(thresholds.size(), false);
  std::vector<uint32_t> n_reads(thresholds.size(), 0);
  MaskProbs(line, target_groups, &mask);
  uint32_t ec_id = 0;
  while (std::getline(probs_file, line)) {
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ','); // first element is the ec id
    uint32_t ref_id = 0;
    while(std::getline(partition, part, ',')) {
      long double abundance = std::stold(part);
      (*assignments)[ec_id][ref_id] = (abundance >= thresholds[ref_id]) && mask[ref_id];
      if ((*assignments)[ec_id][ref_id]) {
	for (uint32_t i = 0; i < assigned_reads[ec_id].size(); ++i) {
	  (*bins)[ref_id].push_back(assigned_reads[ec_id][i] + 1);
	}
      }
      ++ref_id;
    }
    ++ec_id;
  }
  return mask;
}

void BinReads(const std::vector<std::vector<bool>> &assignments, const std::vector<bool> &groups_to_assign, const std::vector<std::vector<uint32_t>> &aligned_reads, std::vector<std::vector<uint32_t>> *assigned_reads) {
  // This function DOES NOT WORK and the AssignProbs func will do what
  // this function is supposed to do.
  uint32_t num_ecs = assignments.size();
  uint32_t n_groups = assignments[0].size();
  uint32_t aln_total = 0;
  for (uint32_t j = 0; j < num_ecs; ++j) {
    uint32_t n_aligned_reads = aligned_reads[j].size();
    aln_total += n_aligned_reads;
    uint32_t bin_id = 0;
    for (uint32_t k = 0; k < n_groups; ++k) {
      if (assignments[j][k] && groups_to_assign[k]) {
  	for (uint32_t h = 0; h < n_aligned_reads; ++h) {
  	  (*assigned_reads)[bin_id].emplace_back(aligned_reads[j][h] + 1);
  	}
  	++bin_id;
      }
    }
  }
}

void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of) {
  for (uint32_t i = 0; i < binned_reads.size(); ++i) {
    of << binned_reads[i] << '\n';
  }
  of.flush();
}
}
