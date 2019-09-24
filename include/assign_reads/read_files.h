#ifndef READ_FILES_H
#define READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

std::unordered_map<long unsigned, std::vector<bool>> read_probs(const double &theta_frac, std::istream &probs_file, std::istream &abundances_file, std::vector<std::string> &ref_names);
void reads_in_ec(std::istream &sam_file, std::istream &ec_file, std::unordered_map<long unsigned, std::vector<std::string>> *reads_in_ec_num);
void read_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::vector<std::string>> *assignments);
void read_groups(const std::vector<std::string> &ref_names, std::istream &groups_file, std::vector<short unsigned> *group_indices);

#endif
  
