#ifndef ASSIGN_READS_CONSTRUCT_THRESHOLDS_H
#define ASSIGN_READS_CONSTRUCT_THRESHOLDS_H

#include <vector>
#include <string>
#include <fstream>

void construct_thresholds(const uint64_t num_ecs, std::istream &abundances_file, std::vector<std::pair<std::string, long double>> *thresholds, long double theta_frac);

#endif
