#ifndef ASSIGN_READS_READ_FILES_H
#define ASSIGN_READS_READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <utility>

void read_probs(const std::vector<std::pair<std::string, long double>> &abundances, std::istream &probs_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *probs);
void read_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *assignments);
void read_groups(const std::vector<std::pair<std::string, long double>> &ref_names, std::istream &groups_file, std::vector<short unsigned> *group_indices);
std::vector<std::pair<std::string, long double>> read_abundances(std::istream &abundances_file);

#endif
