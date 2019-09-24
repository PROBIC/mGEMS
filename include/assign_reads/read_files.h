#ifndef READ_FILES_H
#define READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

std::unordered_map<long unsigned, std::vector<bool>> read_probs(const std::string &probs_path, const std::string &abundances_path, std::vector<std::string> &ref_names, const double &theta_frac);
std::unordered_map<long unsigned, std::vector<std::string>> reads_in_ec(std::istream &sam_file, std::istream &ec_file);
std::unordered_map<long unsigned, std::vector<std::string>> read_assignments(const std::string &assignments_path);
std::vector<short unsigned> read_groups(const std::string &groups_path, const std::vector<std::string> &ref_names);

#endif
  
