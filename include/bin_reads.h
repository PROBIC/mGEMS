#ifndef MGEMS_ASSIGN_READS_H
#define MGEMS_ASSIGN_READS_H

#include <fstream>
#include <vector>
#include <string>

#include "telescope.hpp"

namespace mGEMS {
uint32_t ReadAbundances(std::istream &stream, std::vector<long double> *abundances, std::vector<std::string> *groups);

void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of);
void WriteAssignments(const std::vector<std::vector<bool>> &assignments_mat, const ThemistoAlignment &aln, std::ostream &of);

// mGEMS::BinReads
//   Returns a 2D vector that contains the ids (line numbers in the
//   .fastq files divided by 4) of reads assigned to the groups
//   that were requested.
//   Input:
//     `aln`: Pseudoalignments from Themisto.
//     `abundances`: Relative abundances from mSWEEP.
//     `theta_frac`: Tuning parameter for the thresholds..
//     `single_only`: Only assign reads that are assigned to just a single lineage.
//     `probs_file`: Read probability matrix (.probs file) from mSWEEP.
//     `*target_groups`: Names of the groups that bins will be created for.
//   Output:
//     `*target_groups`: The names will be reordered to match the order of the bins.
//      `*unassigned_bin`: Vector containing the ids of reads that were not assigned to any bin.
//      `out_bins`: Vector containing the bins for the groups given in `*target_groups`.
//      `*assignments_mat`: The read assignment matrix from AssignProbs.
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const std::vector<long double> &abundances, const long double theta_frac, const bool single_only, std::istream &probs_file, std::vector<std::string> *target_groups, std::vector<uint32_t> *unassigned_bin, std::vector<std::vector<bool>> *assignments_mat);
}

#endif
