#ifndef BUILD_SAMPLE_ASSIGN_READS_H
#define BUILD_SAMPLE_ASSIGN_READS_H

#include <vector>

#include "file.hpp"

std::vector<unsigned long> read_assignments(std::istream &assignment_file);
void assign_reads(const std::vector<long unsigned> &assignments, File::Out outfiles[2], File::In infiles[2]);

#endif
