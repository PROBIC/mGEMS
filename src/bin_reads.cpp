#include "assign_reads.h"

#include <sstream>
#include <iterator>
#include <algorithm>

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

void AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, std::vector<std::vector<bool>> *assignments) {
  uint32_t n_groups = thresholds.size();
  std::string line;
  std::getline(probs_file, line); // 1st line is header
  uint32_t ec_id = 0;
  while (std::getline(probs_file, line)) {
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ','); // first element is the ec id

    uint32_t ref_id = 0;
    while(std::getline(partition, part, ',')) {
      long double abundance = std::stold(part);
      (*assignments)[ec_id][ref_id] = (abundance >= thresholds[ref_id]);
      ++ref_id;
    }
    ++ec_id;
  }
}

void BinReads(const std::vector<std::vector<bool>> &assignments, const std::vector<bool> &groups_to_assign, const std::vector<std::vector<uint32_t>> &aligned_reads, std::vector<std::vector<uint32_t>> *assigned_reads) {
  uint32_t num_ecs = assignments.size();
  uint32_t n_groups = assignments[0].size();
  for (uint32_t i = 0; i < num_ecs; ++i) {
    uint32_t group_id = 0;
    for (uint32_t j = 0; j < n_groups; ++j) {
      if (assignments[i][j] && groups_to_assign[j]) {
	for (uint32_t k = 0; k < aligned_reads[i].size(); ++k) {
	  (*assigned_reads)[group_id].emplace_back(aligned_reads[i][k] + 1);
	}
	++group_id;
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
