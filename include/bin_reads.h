#ifndef MGEMS_ASSIGN_READS_H
#define MGEMS_ASSIGN_READS_H

#include <fstream>
#include <vector>
#include <string>

#include "telescope.hpp"

namespace mGEMS {
uint32_t ReadAbundances(std::istream &stream, std::vector<long double> *abundances, std::vector<std::string> *groups);
void ConstructThresholds(const uint32_t num_ecs, const long double theta_frac, const std::vector<long double> &abundances, std::vector<long double> *thresholds);
 std::vector<bool> AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, std::vector<std::string> *target_groups, std::vector<std::vector<bool>> *assignments, const std::vector<std::vector<uint32_t>> &assigned_reads, std::vector<std::vector<uint32_t>> *bins);
void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of);
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::istream &probs_file, std::vector<std::string> *target_groups);
void Extract(const std::vector<std::string> &target_groups, const std::string &outdir, const std::string &strand_1, const std::string &strand_2, std::vector<std::vector<uint32_t>> &bins);
}

#endif
