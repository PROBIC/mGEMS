#ifndef READ_FILES_H
#define READ_FILES_H

#include <string>
#include <unordered_map>
#include <vector>

std::unordered_map<std::vector<bool>, std::vector<std::string>> read_sam(const std::string &sam_path);
std::unordered_map<long unsigned, std::vector<bool>> read_probs(const std::string &probs_path, const std::string &abundances_path, std::vector<std::string> &ref_names);
std::unordered_map<std::vector<bool>, long unsigned> read_ec_ids(const std::string &ec_path, const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec);
  
#endif
  
