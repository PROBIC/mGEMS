#ifndef MGEMS_ASSIGN_READS_H
#define MGEMS_ASSIGN_READS_H

#include <fstream>
#include <vector>
#include <string>

#include "telescope.hpp"

namespace mGEMS {
uint32_t ReadAbundances(std::istream &stream, std::vector<long double> *abundances, std::vector<std::string> *groups);

// Constructs the thresholds according to Equation (7) in the mGEMS
// manuscript. Output will be stored in `thresholds`
// Input:
//   `num_ecs`: How many equivalence classes (unique pseudoalignments) there are in total.
//   `theta_frac`: Tuning parameter for the thresholds. Can be used to loosen/tighten the
//                 rule defined in Equation (7). Should be between 0 and 1 (not validated!).
//   `abundances`: Relative abundances from mSWEEP.
//   `*thresholds`: Holds the output values, i. e. the thresholds. Should be initialized to
//                  the correct size (number of groups in `abundances`) by the caller.
void ConstructThresholds(const uint32_t num_ecs, const long double theta_frac, const std::vector<long double> &abundances, std::vector<long double> *thresholds);

// Performs the actual binning based on the precaculated thresholds.
// Input:
//   `thresholds`: The binning thresholds from CalculateThresholds.
//    `probs_file`: Read probability matrix (.probs file) from mSWEEP.
//    `mask`: Boolean vector defining which groups (value 1) to perform binning on.
//    `*assignments`: `num_ecs` x `n_groups` boolean matrix that contains a 1 if
//                    the ec corresponding to the row was assigned to the column.
//    `aligned_reads`: 2D vector containing the ids of the pseudoaligned reads in
//                     each equivalence class
//    `bins`: 2D output vector containing the ids of the reads binned to each group.
void AssignProbs(const std::vector<long double> &thresholds, std::istream &probs_file, const std::vector<bool> &mask, std::vector<std::vector<bool>> *assignments, const std::vector<std::vector<uint32_t>> &aligned_reads, std::vector<std::vector<uint32_t>> *bins);

void WriteBin(const std::vector<uint32_t> &binned_reads, std::ostream &of);

// Returns a 2D vector that contains the ids (line numbers in the
// .fastq files divided by 4) of reads assigned to the groups
// that were requested.
// Input:
//   `aln`: Pseudoalignments from Themisto.
//   `theta_frac`: Tuning parameter for the thresholds,
//                 should be between 0 and 1 (not validated!).
//   `abundances`: Relative abundances from mSWEEP.
//   `probs_file`: Read probability matrix (.probs file) from mSWEEP.
//   `*target_groups`: Names of the groups that bins will be created for.
// Output:
//   `out_bins`: Vector containing the bins for the groups given in `*target_groups`.
std::vector<std::vector<uint32_t>> Bin(const ThemistoAlignment &aln, const long double theta_frac, const std::vector<long double> &abundances, std::istream &probs_file, std::vector<std::string> *target_groups);

void Extract(const std::vector<std::string> &target_groups, const std::string &outdir, const std::string &strand_1, const std::string &strand_2, std::vector<std::vector<uint32_t>> &bins);
}

#endif
