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
#include "zstr/zstr.hpp"

std::vector<long double> read_abundances(std::istream &abundances_file, std::vector<std::string> &ref_names) {
  std::vector<long double> abundances;
  if (abundances_file.good()) {
    std::string line;
    while (getline(abundances_file, line)) {
      if (line.at(0) != '#') { // skip header lines
	bool abundance = false;
	std::string part;
	std::stringstream partition(line);
	while (getline(partition, part, '\t')) {
	  if (abundance) {
	    abundances.emplace_back(std::stold(part));
	  } else {
	    ref_names.push_back(part);
	    abundance = true;
	  }
	}
      }
    }
  }
  return abundances;
}

void read_probs(const std::vector<long double> &abundances, std::istream &probs_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *ec_to_cluster) {
  uint16_t num_refs = abundances.size();
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
	  std::vector<bool> refs(num_refs, false);
	  ec_to_cluster->at(ec_id).second.resize(num_refs);
	  std::fill(ec_to_cluster->at(ec_id).second.begin(), ec_to_cluster->at(ec_id).second.end(), false);
	  ++ref_id;
	} else {
	  long double abundance = std::stold(part);
	  ec_to_cluster->at(ec_id).second[ref_id - 1] = (abundance >= abundances.at(ref_id - 1));
	  ++ref_id;
	}
      }
    }
  }
}

void read_assignments(std::istream &assignments_file, std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> *assignments) {
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
	  assignments->insert(std::make_pair(current_ec_id, std::make_pair(reads, assign_to)));
	  at_ec_id = false;
	} else {
	  assignments->at(current_ec_id).first.push_back(part);
	}
      }
    }
  }
}

void read_groups(const std::vector<std::string> &ref_names, std::istream &groups_file, std::vector<short unsigned> *group_indices) {
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
  for (size_t i = 0; i < ref_names.size(); ++i) {
    if (groups.find(ref_names.at(i)) != groups.end()) {
      group_indices->push_back(i);
    }
  }
}
