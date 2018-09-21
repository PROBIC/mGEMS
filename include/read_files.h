#ifndef READ_FILES_H
#define READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>

std::unordered_map<long unsigned, std::vector<bool>> read_probs(const std::string &probs_path, const std::string &abundances_path, std::vector<std::string> &ref_names);
std::unordered_map<long unsigned, std::vector<std::string>> reads_in_ec(const std::string &sam_file, const std::string &ec_file);
std::unordered_map<long unsigned, std::vector<std::string>> read_assignments(const std::string &assignments_path);

#endif
  
