#include "assign_reads/construct_thresholds.h"

#include <sstream>
#include <cmath>

void multiply_abundances(std::vector<std::pair<std::string, long double>> *thresholds, long double log_thresh) {
  for (size_t i = 0; i < thresholds->size(); ++i) {
    thresholds->at(i).second = std::exp(std::log(thresholds->at(i).second) + log_thresh);
  }
}

void read_abundances(std::istream &abundances_file, std::vector<std::pair<std::string, long double>> *abundances) {
  if (abundances_file.good()) {
    std::string line;
    while (getline(abundances_file, line)) {
      if (line.at(0) != '#') { // skip header lines
	abundances->emplace_back(std::pair<std::string, long double>());
	bool at_abundance = false;
	std::string part;
	std::stringstream partition(line);
	while (getline(partition, part, '\t')) {
	  if (at_abundance) {
	    abundances->back().second = std::stold(part);
	  } else {
	    abundances->back().first = part;
	    at_abundance = true;
	  }
	}
      }
    }
  }
}

void construct_thresholds(const uint64_t num_ecs, std::istream &abundances_file, std::vector<std::pair<std::string, long double>> *thresholds) {
  read_abundances(abundances_file, thresholds);
  long double log_thresh = std::log1pl(-(long double)thresholds->size()/(long double)num_ecs);
  multiply_abundances(thresholds, log_thresh);
}

void construct_thresholds(const uint64_t num_ecs, std::istream &abundances_file, std::vector<std::pair<std::string, long double>> *thresholds, long double theta_frac) {
  read_abundances(abundances_file, thresholds);
  long double log_thresh = std::log1pl(-(long double)thresholds->size()/(long double)num_ecs);
  log_thresh += theta_frac;
  multiply_abundances(thresholds, log_thresh);
}
