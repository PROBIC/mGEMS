#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <utility>

#include "assign_reads/read_files.h"

void assign_reads(const std::vector<std::pair<std::string, long double>> &thresholds, std::istream &probs_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *reads_in_ec) {
  uint16_t num_refs = thresholds.size();
  if (probs_file.good()) {
    std::string line;
    getline(probs_file, line); // 1st line is header
    uint64_t linenr = 0;
    while (getline(probs_file, line)) {
      std::string part;
      std::stringstream partition(line);
      uint16_t ref_id = 0;
      uint64_t ec_id;
      ++linenr;

      while(getline(partition, part, ',')) {
	if (ref_id == 0) {
	  ec_id = std::stoul(part);
	  reads_in_ec->at(ec_id).second.resize(num_refs);
	  std::fill(reads_in_ec->at(ec_id).second.begin(), reads_in_ec->at(ec_id).second.end(), false);
	  ++ref_id;
	} else {
	  long double abundance = std::stold(part);
	  reads_in_ec->at(ec_id).second[ref_id - 1] = (abundance >= thresholds.at(ref_id - 1).second);
	  ++ref_id;
	}
      }
    }
  }
}

void read_ec_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *reads_in_ec) {
  if (assignments_file.good()) {
    std::string line;
    while(getline(assignments_file, line)) {
      std::string part;
      std::stringstream partition(line);
      bool at_ec_id = true;
      uint64_t current_ec_id;
      while (getline(partition, part, ',')) {	
	if (at_ec_id) {
	  current_ec_id = std::stoul(part);
	  std::vector<std::string> reads;
	  std::vector<bool> assign_to;
	  reads_in_ec->insert(std::make_pair(current_ec_id, std::make_pair(reads, assign_to)));
	  at_ec_id = false;
	} else {
	  reads_in_ec->at(current_ec_id).first.push_back(part);
	}
      }
    }
  }
}

void read_groups_filter(const std::vector<std::pair<std::string, long double>> &thresholds, std::istream &groups_file, std::vector<short unsigned> *group_indices) {
  std::set<std::string> groups;
  if (groups_file.good()) {
    std::string line;
    while(getline(groups_file, line)) {
      std::string part;
      std::stringstream partition(line);
      while (getline(partition, part, '\t')) {
	groups.insert(part);
      }
    }
  }
  for (size_t i = 0; i < thresholds.size(); ++i) {
    if (groups.find(thresholds.at(i).first) != groups.end()) {
      group_indices->push_back(i);
    }
  }
}
