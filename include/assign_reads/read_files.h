#ifndef ASSIGN_READS_READ_FILES_H
#define ASSIGN_READS_READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

void read_probs(const std::vector<double> &abundances, std::istream &probs_file, std::unordered_map<long unsigned, std::vector<bool>> *probs, uint64_t num_ecs);
void reads_in_ec(std::istream &sam_file, const std::string &ec_path, std::unordered_map<long unsigned, std::vector<std::string>> *reads_in_ec_num);
void read_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::vector<std::string>> *assignments, uint64_t &num_ecs);
void read_groups(const std::vector<std::string> &ref_names, std::istream &groups_file, std::vector<short unsigned> *group_indices);
std::vector<double> read_abundances(std::istream &abundances_file, std::vector<std::string> &ref_names);

#endif
