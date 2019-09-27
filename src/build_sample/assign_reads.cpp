#include "build_sample/assign_reads.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>

std::map<std::string, std::set<short unsigned>> read_assignments(std::istream &assignment_file, const unsigned short assignment_id) {
  std::map<std::string, std::set<short unsigned>> assignments;

  if (assignment_file.good()) {
    std::cout << "Reading assignment file" << std::endl;
    std::string read_id; // input file should contain only read names ** as they appear in the fastq files **
    while (getline(assignment_file, read_id)) {
      read_id = "@" + read_id;
      if (assignments.find(read_id) == assignments.end()) {
	std::set<short unsigned> assign;
	assignments.insert(std::pair<std::string, std::set<short unsigned>>(read_id, assign));
      }
      assignments.at(read_id).insert(assignment_id);
    }
  }
  return assignments;
}

void assign_reads(std::unique_ptr<std::ostream> outfiles[1][2], std::istream &strand_1, std::istream &strand_2, const bool gzip_output, const std::map<std::string, std::set<short unsigned>> &assignments) {
  unsigned short K = 1;

  std::string line;
  long unsigned line_nr = 0;
  while (getline(strand_1, line)) {
    if ((line_nr % 4) == 0) {
      std::stringstream parts(line);
      std::string part;
      bool read_name = true;
      while (getline(parts, part, ' ') && read_name) {
  	read_name = false;
  	if (assignments.find(part) != assignments.end()) {
  	  std::string read_id = part;
  	  for (short unsigned val : assignments.at(read_id)) {
  	    *outfiles[val][0] << line << '\n';	  
  	  }
  	  for (short unsigned j = 0; j < 3; ++j) {
  	    getline(strand_1, line);
  	    ++line_nr;
  	    for (short unsigned val : assignments.at(read_id)) {
  	      *outfiles[val][0] << line << '\n';	  
  	    }
  	  }
  	}
      }
    }
    ++line_nr;
  }
  for (size_t i = 0; i < K; ++i) {
    outfiles[i][0]->flush();
  }
  
  line_nr = 0;
  while (getline(strand_2, line)) {
    if ((line_nr % 4) == 0) {
      std::stringstream parts(line);
      std::string part;
      bool read_name = true;
      while (getline(parts, part, ' ') && read_name) {
  	read_name = false;
  	if (assignments.find(part) != assignments.end()) {
  	  std::string read_id = part;
  	  for (short unsigned val : assignments.at(read_id)) {
  	    *outfiles[val][1] << line << '\n';	  
  	  }
  	  for (short unsigned j = 0; j < 3; ++j) {
  	    getline(strand_2, line);
  	    ++line_nr;
  	    for (short unsigned val : assignments.at(read_id)) {
  	      *outfiles[val][1] << line << '\n';	  
  	    }
  	  }
  	}
      }
    }
    ++line_nr;
  }
  for (size_t i = 0; i < K; ++i) {
    outfiles[i][1]->flush();
  }
}
