#include "build_sample/assign_reads_themisto.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>

std::vector<unsigned long> read_assignments(std::istream &assignment_file) {
  std::vector<unsigned long> assignments;

  if (assignment_file.good()) {
    std::cout << "Reading assignment file" << std::endl;
    std::string read_id; // input file should contain only read names ** as they appear in the fastq files **
    while (getline(assignment_file, read_id)) {
	assignments.emplace_back(std::stoul(read_id) + 1);
    }
  }
  std::sort(assignments.rbegin(), assignments.rend());
  return assignments;
}

void process_strand(std::istream &instrand, std::ostream &outstrand, std::vector<long unsigned> assignments) {
  std::string line;
  long unsigned line_nr = 0;
  while (getline(instrand, line)) {
      ++line_nr;
      long unsigned read_id = (line_nr - 1)/4 + 1;
      if(assignments.back() == read_id) {
	  outstrand << line << '\n';
	  for (uint8_t j = 0; j < 3; ++j) {
	      getline(instrand, line);
	      ++line_nr;
	      outstrand << line << '\n';
	  }
	  assignments.pop_back();
      } else {
	  for (uint8_t j = 0; j < 3; ++j) {
	      getline(instrand, line);
	      ++line_nr;
	  }
      }
  }
  outstrand.flush();
}

void assign_reads(const std::vector<long unsigned> &assignments, File::Out outfiles[2], File::In infiles[2]) {
    process_strand(infiles[0].stream(), outfiles[0].stream(), assignments);
    process_strand(infiles[1].stream(), outfiles[1].stream(), assignments);
}
