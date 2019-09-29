#ifndef ASSIGN_READS_WRITE_FILES_H
#define ASSIGN_READS_WRITE_FILES_H

#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>

void write_reads(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<short unsigned> &group_indices,
		 const std::vector<std::unique_ptr<std::ostream>> &outfiles);
void write_ecs(const std::unordered_map<long unsigned,
	       std::vector<std::string>> &reads_in_ec,
	       std::unique_ptr<std::ostream> &of);
#endif
