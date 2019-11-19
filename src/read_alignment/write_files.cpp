#include "read_alignment/write_files.h"

void write_ecs(const std::unordered_map<long unsigned, std::vector<std::string>> &reads_in_ec, File::Out &of) {
  if (of.stream().good()) {
    for (auto kv : reads_in_ec) {
      of << kv.first;
      for (auto read : kv.second) {
	of << ',' << read;
      }
      of << '\n';
    }
  }
  of.close();
}
