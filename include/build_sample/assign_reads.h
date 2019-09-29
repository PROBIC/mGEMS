#ifndef BUILD_SAMPLE_ASSIGN_READS_H
#define BUILD_SAMPLE_ASSIGN_READS_H

#include <map>
#include <string>
#include <set>
#include <fstream>
#include <memory>

std::map<std::string, std::set<short unsigned>> read_assignments(std::istream &assignment_file, const unsigned short assignment_id);
void assign_reads(std::unique_ptr<std::ostream> outfiles[1][2], std::unique_ptr<std::istream> infiles[2], const bool gzip_output, const std::map<std::string, std::set<short unsigned>> &assignments);

#endif
