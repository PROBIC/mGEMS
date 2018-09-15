#ifndef WRITE_FILES_H
#define WRITE_FILES_H

#include <map>
#include <unordered_map>
#include <vector>
#include <string>

void write_reads(const std::map<std::vector<bool>, long unsigned> &ec_to_id, const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec, std::unordered_map<long unsigned, std::vector<bool>> probs, std::vector<std::string> ref_names, const std::string &outfile);
#endif
