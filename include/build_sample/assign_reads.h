#ifndef ASSIGN_READS_H
#define ASSIGN_READS_H

#include <map>
#include <string>
#include <set>

std::map<std::string, std::set<short unsigned>> read_assignments(const std::string &assignment_path, const unsigned short assignment_id);
void assign_reads(const std::string &outfile, const std::string &strand1, const std::string &strand2, const bool gzip_output, const std::map<std::string, std::set<short unsigned>> &assignments);

#endif
