#ifndef ASSIGN_READS_READ_FILES_H
#define ASSIGN_READS_READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <utility>

void assign_reads(const std::vector<std::pair<std::string, long double>> &thresholds, std::istream &probs_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *reads_in_ec);
void read_ec_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *reads_in_ec);
void read_groups_filter(const std::vector<std::pair<std::string, long double>> &thresholds, std::istream &groups_file, std::vector<short unsigned> *group_indices);

#endif
