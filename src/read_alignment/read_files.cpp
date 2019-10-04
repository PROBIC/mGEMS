#include "read_alignment/read_files.h"

#include <sstream>

#include "zstr/zstr.hpp"

void read_header(const std::string &header_line, std::unordered_map<std::string, long unsigned> &ref_to_id, long unsigned &ref_id) {
  if (header_line.at(1) == 'S') {
    std::stringstream parts(header_line);
    std::string part;
    
    short firstel = 0;
    while (getline(parts, part, '\t')) {
      if (firstel == 1) {
	ref_to_id[part.substr(3)] = ref_id;
	++ref_id;
      }
      ++firstel;
    }
  }
}

void read_alignment(const std::string &alignment_line, const std::unordered_map<std::string, long unsigned> &ref_to_id, std::unordered_map<std::string, std::vector<bool>> &read_to_ec) {
  std::stringstream parts(alignment_line);
  std::string part;
  unsigned n_refs = ref_to_id.size();
  
  short firstel = 0;
  std::string read_id;
  while (getline(parts, part, '\t') && firstel < 3) { // we only want the read (0) and reference (2)
    if (firstel == 0) {
      read_id = part;
      if (read_to_ec.find(part) == read_to_ec.end()) {
	std::vector<bool> ec_config(n_refs, 0);
	read_to_ec[part] = ec_config;
      }
    } else if (firstel == 2) {
      if (part != "*") { // star signifies no alignment
	read_to_ec[read_id][ref_to_id.at(part)] = 1;
      }
    }
    ++firstel;
  }
}

std::unordered_map<std::vector<bool>, std::vector<std::string>> read_sam(std::istream &sam_file) {
  std::unordered_map<std::string, long unsigned> ref_to_id;
  std::unordered_map<std::string, std::vector<bool>> read_to_ec;
  std::unordered_map<std::vector<bool>, std::vector<std::string>> reads_in_ec;
  if (sam_file.good()) {
    std::string line;
    long unsigned ref_id = 0;
    while (getline(sam_file, line)) {
      if (line.at(0) == '@') {
	read_header(line, ref_to_id, ref_id);
      } else {
	read_alignment(line, ref_to_id, read_to_ec);
      }
    }
  }

  // Assign reads to equivalence classes.
  for (auto kv : read_to_ec) {
    if (reads_in_ec.find(kv.second) == reads_in_ec.end()) {
      std::vector<std::string> reads;
      reads_in_ec[kv.second] = reads;
    }
    reads_in_ec[kv.second].emplace_back(kv.first);
  }

  return reads_in_ec;
}

std::unordered_map<std::vector<bool>, long unsigned> read_ec_ids(std::istream &ec_file, const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec) {
  std::unordered_map<std::vector<bool>, long unsigned> ec_to_id;
  unsigned n_refs = reads_in_ec.begin()->first.size();
  if (ec_file.good()) {
    std::string line;
    while (getline(ec_file, line)) {
      std::string part;
      std::stringstream partition(line);
      bool firstel = true;
      long unsigned key = 0;
      std::vector<bool> config(n_refs, 0);
      while (getline(partition, part, '\t')) {
	if (firstel) {
	  key = std::stoi(part);
	  firstel = false;
	} else {
	  std::string one;
	  std::stringstream ones(part);
	  while (getline(ones, one, ',')) {
	    unsigned makeone = std::stoi(one);
	    config[makeone] = 1;
	  }
	}
      }
      if (reads_in_ec.find(config) != reads_in_ec.end()) {
        ec_to_id[config] = key;
      }
    }
  }
  return ec_to_id;
}

void reads_in_ec(std::istream &sam_file, const std::string &ec_path, std::unordered_map<long unsigned, std::vector<std::string>> *reads_in_ec_num) {
  const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec = read_sam(sam_file);
  zstr::ifstream ec_file(ec_path);
  const std::unordered_map<std::vector<bool>, long unsigned> &ec_to_id = read_ec_ids(ec_file, reads_in_ec);
  reads_in_ec_num->reserve(ec_to_id.bucket_count());
  for (auto kv : reads_in_ec) {
    if (ec_to_id.find(kv.first) != ec_to_id.end()) {
      reads_in_ec_num->insert(make_pair(ec_to_id.at(kv.first), kv.second));
    }
  }
}
