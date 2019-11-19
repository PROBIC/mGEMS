#ifndef BUILD_SAMPLE_ASSIGN_READS_H
#define BUILD_SAMPLE_ASSIGN_READS_H

#include <string>
#include <set>
#include <fstream>
#include <memory>

#include "file.hpp"

std::set<std::string> read_assignments(std::istream &assignment_file);
void assign_reads(const std::set<std::string> &assignments, File::Out outfiles[2], File::In infiles[2]);

#endif
