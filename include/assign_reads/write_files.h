#ifndef ASSIGN_READS_WRITE_FILES_H
#define ASSIGN_READS_WRITE_FILES_H

#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>

void write_reads(const std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> &reads_in_ec,
		 const std::vector<short unsigned> &group_indices,
		 const std::vector<std::unique_ptr<std::ostream>> &outfiles);
#endif
