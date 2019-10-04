#include "assign_reads/write_files.h"

#include <set>
#include <iostream>

void write_reads(const std::unordered_map<long unsigned, std::pair<std::vector<std::string>, std::vector<bool>>> &reads_in_ec,
		 const std::vector<short unsigned> &group_indices,
		 const std::vector<std::unique_ptr<std::ostream>> &ofs) {
  for (auto kv : reads_in_ec) {
    for (auto i : group_indices) {
      if (kv.second.second.at(i)) {
	for (auto read : kv.second.first) {
	  *ofs[i] << read << '\n';
	}
      }
    }
  }

  for (auto i : group_indices) {
    ofs[i]->flush();
  }
}
