#include "bin_reads.h"

#include <sstream>
#include <iterator>
#include <algorithm>
#include <utility>
#include <cmath>

namespace mGEMS {
void FilterTargetGroups(const std::vector<std::string> &group_names, const std::vector<double> &abundances, const double min_abundance, std::vector<std::string> *target_groups) {
  uint32_t n_groups = group_names.size();
  for (uint32_t i = 0; i < n_groups; ++i) {
    if (abundances[i] < min_abundance && std::find(target_groups->begin(), target_groups->end(), group_names[i]) != target_groups->end()) {
      target_groups->erase(std::find(target_groups->begin(), target_groups->end(), group_names[i]));
    }
  }
}

void ReadAbundances(std::istream &stream, std::vector<double> *abundances, std::vector<std::string> *groups) {
  std::string line;
  while(std::getline(stream, line)) {
    if (line.at(0) != '#') {
      std::string part;
      std::stringstream parts(line);
      std::getline(parts, part, '\t');
      groups->emplace_back(part);
      std::getline(parts, part, '\t');
      abundances->emplace_back(std::stold(part));
    }
  }
}

void ConstructThresholds(const uint32_t num_ecs, const long double theta_frac, const std::vector<double> &abundances, std::vector<long double> *thresholds, bool logscale = false) {
  // Constructs the thresholds according to Equation (7) in the mGEMS
  // manuscript. Output will be stored in `thresholds`
  // Input:
  //   `num_ecs`: How many equivalence classes (unique pseudoalignments) there are in total.
  //   `theta_frac`: Tuning parameter for the thresholds. Can be used to loosen/tighten the
  //                 rule defined in Equation (7).
  //   `abundances`: Relative abundances from mSWEEP.
  // Output:
  //   `*thresholds`: Holds the output values, i. e. the thresholds. Should be initialized to
  //                  the correct size (number of groups in `abundances`) by the caller.
  long double log_thresh = std::log1pl(-(long double)abundances.size()/(long double)num_ecs);
  log_thresh += std::log(theta_frac);
  uint32_t n_refs = abundances.size();
  for (uint32_t i = 0; i < n_refs; ++i) {
    (*thresholds)[i] = std::log(abundances[i]) + log_thresh;
    if (!logscale)
      (*thresholds)[i] = std::exp((*thresholds)[i]);
  }
}


void ReadHeader(const std::string &header_line, const size_t n_groups, std::vector<std::string> *group_names) {
  std::string part;
  std::stringstream partition(header_line);
  std::getline(partition, part, ',');
  for (size_t i = 0; i < n_groups; ++i) {
    std::getline(partition, part, ',');
    group_names->emplace_back(std::move(part));
  }
}

void MaskProbs(const std::vector<std::string> &all_group_names, std::vector<std::string> *target_groups, std::vector<bool> *mask) {
  std::vector<std::string> ordered_targets(target_groups->size());
  uint32_t ref_id = 0;
  uint32_t target_id = 0;
  for (size_t i = 0; i < all_group_names.size(); ++i) {
    if (std::find(target_groups->begin(), target_groups->end(), all_group_names[i]) != target_groups->end()) {
      (*mask)[ref_id] = true;
      ordered_targets[target_id] = all_group_names[i];
      ++target_id;
    }
    ++ref_id;
  }
  (*target_groups) = ordered_targets;
}

bool EvaluateAssignment(const double abundance, const double threshold, uint32_t *n_assignments, bool *any_assigned) {
      bool assigned = abundance >= threshold;
      (*n_assignments) += assigned;
      (*any_assigned) = (*any_assigned) || assigned;
      return assigned;
}

void InsertAssigned(const size_t ec_id, const bool single_only, const size_t n_assignments, const telescope::ThemistoAlignment &alignment, const std::vector<bool> &mask, std::vector<std::vector<bool>> *assignments, std::vector<std::vector<uint32_t>> *bins, std::vector<uint32_t> *unassigned_bin) {
  if (single_only && n_assignments == 1) {
    for (uint32_t j = 0; j < (*assignments)[ec_id].size(); ++j) {
      if ((*assignments)[ec_id][j] && mask[j]) {
	for (uint32_t i = 0; i < alignment.reads_assigned_to_ec(ec_id).size(); ++i) {
	  (*bins)[j].push_back(alignment.reads_assigned_to_ec(ec_id)[i] + 1);
	}
      }
    }
  } else if (!single_only && n_assignments > 0) {
    for (uint32_t j = 0; j < (*assignments)[ec_id].size(); ++j) {
      if ((*assignments)[ec_id][j] && mask[j]) {
	for (uint32_t i = 0; i < alignment.reads_assigned_to_ec(ec_id).size(); ++i) {
	  (*bins)[j].push_back(alignment.reads_assigned_to_ec(ec_id)[i] + 1);
	}
      }
    }
  } else if (n_assignments == 0) {
    // Send reads to aligned but unassigned bin.
    for (uint32_t i = 0; i < alignment.reads_assigned_to_ec(ec_id).size(); ++i) {
      unassigned_bin->push_back(alignment.reads_assigned_to_ec(ec_id)[i] + 1);
    }
  }
}

  void InsertAssigned(const size_t ec_id, const bool single_only, const size_t n_assignments, const std::vector<std::vector<uint32_t>> &reads_to_ecs, const std::vector<bool> &mask, std::vector<std::vector<bool>> *assignments, std::vector<std::vector<uint32_t>> *bins, std::vector<uint32_t> *unassigned_bin) {
  if (single_only && n_assignments == 1) {
    for (uint32_t j = 0; j < (*assignments)[ec_id].size(); ++j) {
      if ((*assignments)[ec_id][j] && mask[j]) {
	for (uint32_t i = 0; i < reads_to_ecs[ec_id].size(); ++i) {
	  (*bins)[j].push_back(reads_to_ecs[ec_id][i] + 1);
	}
      }
    }
  } else if (!single_only && n_assignments > 0) {
    for (uint32_t j = 0; j < (*assignments)[ec_id].size(); ++j) {
      if ((*assignments)[ec_id][j] && mask[j]) {
	for (uint32_t i = 0; i < reads_to_ecs[ec_id].size(); ++i) {
	  (*bins)[j].push_back(reads_to_ecs[ec_id][i] + 1);
	}
      }
    }
  } else if (n_assignments == 0) {
    // Send reads to aligned but unassigned bin.
    for (uint32_t i = 0; i < reads_to_ecs[ec_id].size(); ++i) {
      unassigned_bin->push_back(reads_to_ecs[ec_id][i] + 1);
    }
  }
}

void AssignProbsMatrix(const std::vector<long double> &thresholds, const seamat::Matrix<double> &probs_mat, const std::vector<bool> &mask, const std::vector<std::vector<uint32_t>> &reads_to_ecs, const bool single_only, std::vector<std::vector<bool>> *assignments, std::vector<std::vector<uint32_t>> *bins, std::vector<uint32_t> *unassigned_bin) {
  // Performs the actual binning based on the precaculated thresholds.
  // Input:
  //   `thresholds`: The binning thresholds from CalculateThresholds.
  //    `probs_mat`: Read probability matrix from mSWEEP.
  //    `mask`: Boolean vector defining which groups (value 1) to perform binning on.
  //    `aligned_reads`: 2D vector containing the ids of the pseudoaligned reads in
  //                     each equivalence class
  //    `single_only`: Only assign reads that are assigned to just a single lineage.
  // Output:
  //    `*assignments`: `num_ecs` x `n_groups` boolean matrix that contains a 1 if
  //                    the ec corresponding to the row was assigned to the column.
  //    `*bins`: 2D output vector containing the ids of the reads binned to each group.
  //    `*unassigned_bin`: Vector containing the ids of reads that were not assigned to any bin.
  std::vector<uint32_t> n_reads(thresholds.size(), 0);

  size_t n_alignments = probs_mat.get_cols();
  size_t n_groups = probs_mat.get_rows();

  for (size_t i = 0; i < n_alignments; ++i) {
    bool any_assigned = false;
    uint32_t n_assignments = 0;
    for (size_t j = 0; j < n_groups; ++j) {
      (*assignments)[i][j] = EvaluateAssignment(probs_mat(j, i), thresholds[j], &n_assignments, &any_assigned);
    }
    InsertAssigned(i, single_only, n_assignments, reads_to_ecs, mask, assignments, bins, unassigned_bin);
  }
}

void AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, const std::vector<bool> &mask, const telescope::ThemistoAlignment &alignment, const bool single_only, std::vector<std::vector<bool>> *assignments, std::vector<std::vector<uint32_t>> *bins, std::vector<uint32_t> *unassigned_bin) {
  // Performs the actual binning based on the precaculated thresholds.
  // Input:
  //   `thresholds`: The binning thresholds from CalculateThresholds.
  //    `probs_file`: Read probability matrix (.probs file) from mSWEEP.
  //    `mask`: Boolean vector defining which groups (value 1) to perform binning on.
  //    `aligned_reads`: 2D vector containing the ids of the pseudoaligned reads in
  //                     each equivalence class
  //    `single_only`: Only assign reads that are assigned to just a single lineage.
  // Output:
  //    `*assignments`: `num_ecs` x `n_groups` boolean matrix that contains a 1 if
  //                    the ec corresponding to the row was assigned to the column.
  //    `*bins`: 2D output vector containing the ids of the reads binned to each group.
  //    `*unassigned_bin`: Vector containing the ids of reads that were not assigned to any bin.
  std::vector<uint32_t> n_reads(thresholds.size(), 0);
  uint32_t ec_id = 0;
  std::string line;
  while (std::getline(probs_file, line) && !line.empty()) {
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ','); // first element is the ec id
    uint32_t ref_id = 0;
    bool any_assigned = false;
    uint32_t n_assignments = 0;
    while(std::getline(partition, part, ',')) {
      (*assignments)[ec_id][ref_id] = EvaluateAssignment(std::stold(part), thresholds[ref_id], &n_assignments, &any_assigned);
      ++ref_id;
    }
    InsertAssigned(ec_id, single_only, n_assignments, alignment, mask, assignments, bins, unassigned_bin);
    ++ec_id;
  }
}

void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of) {
  for (uint32_t i = 0; i < binned_reads.size(); ++i) {
    of << binned_reads[i] << '\n';
  }
  of.flush();
}

void WriteAssignments(const std::vector<std::vector<bool>> &assignments_mat, const telescope::ThemistoAlignment &aln, std::ostream &of) {
  uint32_t n_groups = assignments_mat[0].size();
  uint32_t n_ecs = assignments_mat.size();
  for (uint32_t i = 0; i < n_ecs; ++i) {
    for (uint32_t j = 0; j < aln.reads_assigned_to_ec(i).size(); ++j) {
      of << aln.reads_assigned_to_ec(i)[j] + 1<< '\t';
      for (uint32_t k = 0; k < n_groups; ++k) {
	of << assignments_mat[i][k];
	if (k < n_groups - 1) {
	  of << '\t';
	} else {
	  if (i < n_ecs - 1) {
	    of << '\n';
	  }
	}
      }
    }
  }
  of << std::endl;
}

std::vector<std::vector<uint32_t>> BuildOutBins(const size_t n_groups, std::vector<std::vector<uint32_t>> &read_bins, std::vector<uint32_t> *unassigned_bin) {
  std::vector<std::vector<uint32_t>> out_bins;
  for (uint32_t i = 0; i < n_groups; ++i) {
    if (read_bins[i].size() > 0) {
      out_bins.push_back(std::move(read_bins[i]));
      std::sort(out_bins.back().begin(), out_bins.back().end());
    }
  }
  return out_bins;
}

std::vector<std::vector<uint32_t>> BinFromMatrix(const std::vector<std::vector<uint32_t>> &reads_to_ecs, const std::vector<double> &abundances, const seamat::Matrix<double> &probs_mat, const std::vector<std::string> &all_group_names, std::vector<std::string> *target_groups) {
  uint32_t num_ecs = probs_mat.get_cols(); // Probs mat is of size n_groups x num_ecs
  uint32_t n_groups = abundances.size();
  std::vector<long double> thresholds(n_groups);
  long double theta_frac = 1.0;
  ConstructThresholds(num_ecs, theta_frac, abundances, &thresholds, true);

  std::vector<std::vector<uint32_t>> read_bins(n_groups, std::vector<uint32_t>());

  // Build a binary vector indicating which groups to output
  std::vector<bool> mask(n_groups, false);
  MaskProbs(all_group_names, target_groups, &mask);

  bool single_only = false;
  std::vector<std::vector<bool>> assignments_mat(num_ecs, std::vector<bool>(n_groups, false));
  std::vector<uint32_t> unassigned_bin;
  AssignProbsMatrix(thresholds, probs_mat, mask, reads_to_ecs, single_only, &assignments_mat, &read_bins, &unassigned_bin);

  const std::vector<std::vector<uint32_t>> &out_bins = BuildOutBins(n_groups, read_bins, &unassigned_bin);

  return out_bins;
}

std::vector<std::vector<uint32_t>> Bin(const telescope::ThemistoAlignment &aln, const std::vector<double> &abundances, const long double theta_frac, const bool single_only, std::istream &probs_file, std::vector<std::string> *target_groups, std::vector<uint32_t> *unassigned_bin, std::vector<std::vector<bool>> *assignments_mat) {
  uint32_t num_ecs = aln.n_ecs();
  uint32_t n_groups = abundances.size();
  std::vector<long double> thresholds(n_groups);
  ConstructThresholds(num_ecs, theta_frac, abundances, &thresholds);

  std::vector<std::vector<uint32_t>> read_bins(n_groups, std::vector<uint32_t>());

  // Read in all group names from the probs file
  std::string probs_header;
  std::getline(probs_file, probs_header); // 1st line is header
  std::vector<std::string> all_group_names;
  ReadHeader(probs_header, n_groups, &all_group_names);

  // Build a binary vector indicating which groups to output
  std::vector<bool> mask(n_groups, false);
  MaskProbs(all_group_names, target_groups, &mask);

  AssignProbs(thresholds, probs_file, mask, aln, single_only, assignments_mat, &read_bins, unassigned_bin);

  const std::vector<std::vector<uint32_t>> &out_bins = BuildOutBins(n_groups, read_bins, unassigned_bin);
  std::sort(unassigned_bin->begin(), unassigned_bin->end());

  return out_bins;
}
}
