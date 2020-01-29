#ifndef READ_SAM_READ_FILES_H
#define READ_SAM_READ_FILES_H

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "file.hpp"

void reads_in_ec(std::istream &sam_file, File::In &ec_file, std::unordered_map<long unsigned, std::vector<std::string>> *reads_in_ec_num, bool themisto, uint16_t n_refs = 0);

#endif
