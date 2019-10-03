#ifndef BUILD_SAMPLE_ASSIGN_READS_H
#define BUILD_SAMPLE_ASSIGN_READS_H

#include <string>
#include <set>
#include <fstream>
#include <memory>

std::set<std::string> read_assignments(std::istream &assignment_file);
void assign_reads(const std::set<std::string> &assignments, std::unique_ptr<std::ostream> outfiles[2], std::unique_ptr<std::istream> infiles[2]);

#endif
