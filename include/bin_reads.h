#ifndef MGEMS_ASSIGN_READS_H
#define MGEMS_ASSIGN_READS_H

#include <fstream>
#include <vector>
#include <string>

#include "telescope.hpp"

#include "file.hpp"

namespace mGEMS {
uint32_t ReadAbundances(std::istream &stream, std::vector<long double> *abundances, std::vector<std::string> *groups);
void ConstructThresholds(const uint32_t num_ecs, const long double theta_frac, const std::vector<long double> &abundances, std::vector<long double> *thresholds);
std::vector<bool> AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, std::vector<std::string> *target_groups, std::vector<std::vector<bool>> *assignments);
void BinReads(const std::vector<std::vector<bool>> &assignments, const std::vector<bool> &groups_to_assign, const std::vector<std::vector<uint32_t>> &aligned_reads, std::vector<std::vector<uint32_t>> *assigned_reads);
void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of);
}

#endif
