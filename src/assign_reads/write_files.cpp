#include "assign_reads/write_files.h"

#include <set>
#include <iostream>

void write_reads(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec,
		 const std::unordered_map<long unsigned, std::vector<bool>> &probs,
		 const std::vector<short unsigned> &group_indices,
		 const std::vector<std::unique_ptr<std::ostream>> &ofs) {
  for (auto kv : reads_in_ec) {
    bool any_probs = (probs.find(kv.first) != probs.end());
    if (any_probs) {
      for (auto i : group_indices) {
	if (probs.at(kv.first)[i]) {
	  for (auto read : kv.second) {
	    *ofs[i] << read << '\n';
	  }
	}
      }
    }
  }

  for (auto i : group_indices) {
    ofs[i]->flush();
  }
}

void write_ecs(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec, std::unique_ptr<std::ostream> &of) {
  if (of->good()) {
    for (auto kv : reads_in_ec) {
      *of << kv.first;
      for (auto read : kv.second) {
	*of << ',' << read;
      }
      *of << '\n';
    }
  }
  of->flush();
}
