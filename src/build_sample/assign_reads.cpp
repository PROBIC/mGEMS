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

void process_strand(const std::map<std::string, std::set<short unsigned>> &assignments, std::istream &instrand, std::ostream &outstrand) {
  std::string line;
  long unsigned line_nr = 0;
  while (getline(instrand, line)) {
    if ((line_nr % 4) == 0) {
      std::stringstream parts(line);
      std::string part;
      bool read_name = true;
      while (getline(parts, part, ' ') && read_name) {
  	read_name = false;
  	if (assignments.find(part) != assignments.end()) {
  	  std::string read_id = part;
  	  for (short unsigned val : assignments.at(read_id)) {
  	    outstrand << line << '\n';	  
  	  }
  	  for (short unsigned j = 0; j < 3; ++j) {
  	    getline(instrand, line);
  	    ++line_nr;
  	    for (short unsigned val : assignments.at(read_id)) {
  	      outstrand << line << '\n';	  
  	    }
  	  }
  	}
      }
    }
    ++line_nr;
  }
    outstrand.flush();
}

void assign_reads(const std::map<std::string, std::set<short unsigned>> &assignments, std::unique_ptr<std::ostream> outfiles[2], std::unique_ptr<std::istream> infiles[2]) {
  process_strand(assignments, *infiles[0], *outfiles[0]);
  process_strand(assignments, *infiles[1], *outfiles[1]);
}
