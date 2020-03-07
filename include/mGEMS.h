#ifndef MGEMS_MGEMS_H
#define MGEMS_MGEMS_H

#include <vector>
#include <string>
#include <fstream>

#include "telescope.hpp"

namespace mGEMS {
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::vector<std::string> &group_names, std::istream &probs_file, std::vector<std::string> *target_groups);
void Extract(const std::vector<std::string> &target_groups, const std::string &outdir, const std::string &strand_1, const std::string &strand_2, std::vector<std::vector<uint32_t>> &bins);
}

#endif
