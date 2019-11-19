#ifndef READ_SAM_WRITE_FILES_H
#define READ_SAM_WRITE_FILES_H

#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>

#include "file.hpp"

void write_ecs(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec, File::Out &of);

#endif
