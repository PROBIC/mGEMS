#ifndef READ_SAM_READ_FILES_H
#define READ_SAM_READ_FILES_H

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

void reads_in_ec(std::istream &sam_file, const std::string &ec_path, std::unordered_map<long unsigned, std::vector<std::string>> *reads_in_ec_num);

#endif
