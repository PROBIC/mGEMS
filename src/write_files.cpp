#include <fstream>
#include <set>

#include "write_files.h"

void write_ec_to_id(const std::map<std::vector<bool>, long unsigned> &ec_to_id, const std::unordered_map<std::vector<bool>, std::vector<std::string>> &reads_in_ec, const std::string &outfile) {
  std::ofstream of(outfile);

  for (auto kv : reads_in_ec) {
    if (ec_to_id.find(kv.first) == ec_to_id.end()) {
      unsigned abc = 0;
      for (auto val : kv.first) {
	abc += val;
      }
    } else {
      of << ec_to_id.at(kv.first) << '\t';
    }
    for (long unsigned i = 0; i < kv.second.size(); ++i) {
      of << kv.second[i];
      if (i != kv.second.size()) {
	of << '\t';
      }
    }
    of << '\n';
  }
  of.close();
}
