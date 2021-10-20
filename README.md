# mGEMS

Bacterial sequencing data binning on strain-level based on probabilistic taxonomic classification.

More about mGEMS in the preprint [Genomic epidemiology with mixed
samples](https://www.biorxiv.org/content/10.1101/2020.04.03.021501v2)
in bioRxiv (not peer reviewed).

## Installation
### Dependencies
To run the binning + assembly pipeline, you will need a program that
does pseudoalignment and another program that estimates an assignment
probability matrix for the reads to the alignment targets.

We recommend to use [Themisto](https://github.com/algbio/themisto)
(v0.1.1 or newer) for pseudoalignment and
[mSWEEP](https://github.com/probic/mSWEEP) (v1.3.2 or newer) for
estimating the probability matrix. For assembling the bins output by
mGEMS, we recommend [shovill](https://github.com/tseemann/shovill) for
typical use-cases but metagenomic assemblers like
[MEGAHIT](https://github.com/voutcn/megahit) may perform better when
the differences between the bins are especially small (see
Supplementary Figure 2 of the mGEMS preprint).

### mGEMS binaries
mGEMS can be obtained in the form of a precompiled binary
* [Linux 64-bit binary](https://github.com/PROBIC/mGEMS/releases/download/v1.0.0/mGEMS_linux-v1.1.0.tar.gz)
* [macOS 64-bit binary](https://github.com/PROBIC/mGEMS/releases/download/v1.0.0/mGEMS_macOS-v1.1.0.tar.gz)

or by following the instructions below for compiling mGEMS from source.

### Compiling from source
#### Requirements
- C++11 compliant compiler.
- cmake

#### Compilation
Clone the repository
```
git clone https://github.com/PROBIC/mGEMS.git
```
enter the directory and run
```
mkdir build
cd build
cmake ..
make
```
This will compile the mGEMS executable in the build/bin/ directory.

## Usage
### mGEMS
The mGEMS executable provides three commands: mGEMS, mGEMS bin, and
mGEMS extract. The first command (mGEMS) is shorthand for running both
mGEMS bin and mGEMS extract, which bin the reads in the input
pseudoalignment (mGEMS bin) and extract the binned reads from the
original mixed samples (mGEMS extract).

### Tutorial — E. coli ST131 sublineages
A tutorial for reproducing the *E. coli* ST131 sublineage phylogenetic
tree presented in Mäklin et al. 2020 using mGEMS is available in the
[docs folder of this repository](docs/TUTORIAL.md).

### Quickstart — full pipeline
Build a [Themisto](https://github.com/algbio/themisto) index to
align against.
```
mkdir themisto_index
mkdir themisto_index/tmp
build_index --k 31 --input-file example.fasta --auto-colors --index-dir themisto_index --temp-dir themisto_index/tmp
```

Align paired-end reads 'reads_1.fastq.gz' and 'reads_2.fastq.gz' with Themisto (note the **--sort-output** flag must be used!)
```
pseudoalign --index-dir themisto_index --query-file reads_1.fastq.gz --outfile pseudoalignments_1.aln --rc --temp-dir themisto_index/tmp --n-threads 16 --mem-megas 8192 --sort-output --gzip-output
pseudoalign --index-dir themisto_index --query-file reads_2.fastq.gz --outfile pseudoalignments_2.aln --rc --temp-dir themisto_index/tmp --n-threads 16 --mem-megas 8192 --sort-output --gzip-output
```

Estimate the relative abundances with mSWEEP (reference_grouping.txt
should contain the groups the sequences in 'example.fasta' are
assigned to. See the [mSWEEP](https://github.com/probic/mSWEEP) usage instructions for details).
```
mSWEEP --themisto-1 pseudoalignments_1.aln.gz --themisto-2 pseudoalignments_2.aln.gz -o mSWEEP -i reference_grouping.txt --write-probs
```

Bin the reads and write all bins to the 'mGEMS-out' folder
```
mkdir mGEMS-out
mGEMS -r reads_1.fastq.gz,reads_2.fastq.gz --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.csv -a mSWEEP_abundances.txt --index themisto_index
```
This will write the binned paired-end reads for *all groups* in the
mSWEEP_abundances.txt file in the mGEMS-out folder (compressed with
zlib).

### Advanced use
You can also extract the read-to-group assignments table that mGEMS
uses internally by adding the `--write-assignment-table` toggle to the
call to `mGEMS` or `mGEMS bin`:
```
mGEMS --groups group-3,group-4 -r reads_1.fastq.gz,reads_2.fastq.gz --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.csv -a mSWEEP_abundances.txt --index themisto_index --write-assignment-table
```

... or bin and write only the reads that are assigned to "group-3" or
"group-4" by adding the '--groups group-3,group-4' flag
```
mGEMS --groups group-3,group-4 -r reads_1.fastq.gz,reads_2.fastq.gz --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.csv -a mSWEEP_abundances.txt --index themisto_index
```

... write the reads that pseudoaligned to a reference sequence but were not assigned to any group by adding the `--write-unassigned` flag:
```
mGEMS --groups group-3,group-4 -r reads_1.fastq.gz,reads_2.fastq.gz --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.csv -a mSWEEP_abundances.txt --index themisto_index --write-unassigned
```

Alternatively, find and write only the read bins for "group-3",
"group-4", and the reads that pseudoaligned but were not assigned to
any group; skipping extracting the reads
```
mGEMS bin --groups group-3,group-4 --themisto-alns pseudoalignments_1.aln.gz,pseudoalignments_2.aln.gz -o mGEMS-out --probs mSWEEP_probs.csv -a mSWEEP_abundances.txt --index themisto_index --write-unassigned
```

... and extract the reads when feeling like it
```
mGEMS extract --bins mGEMS-out/group-3.bin,mGEMS-out/group-4.bin,mGEMS-out/unassigned_reads.bin -r
reads_1.fastq.gz,reads_2.fastq.gz -o mGEMS-out
```

### Accepted input flags
mGEMS accepts the following input flags
```
	-r                       Comma-separated list of input read(s).
	--themisto-alns          Comma-separated list of pseudoalignment file(s) 
	                         for the reads from themisto.
	-o                       Output directory (must exist before running!).
	--probs                  Comma-separated Posterior probability matrix (output from mSWEEP with
	                         the --write-probs flag).
	-a                       Relative abundance estimates from mSWEEP (tab-separated, 1st
	                         column has the group names and 2nd column the estimates).
	--index                  Themisto pseudoalignment index directory.
	--groups                 (Optional) Which groups to extract from the input reads.
	--min-abundance          (Optional) Extract only groups that have a relative abundance higher than this value.
	--compress               (Optional) Toggle compressing the output files (default: compress)
	--write-unassigned       (Optional) Extract reads that pseudoaligned to a reference sequence but were not assigned to any group. (default: off)
	--write-assignment-table (Optional) Write the read to group assignments table to `reads_to_groups.tsv` in the output directory. (default: off).
	--unique-only            (Optional) Write only the reads that are assigned to a single group.
```

## License
The source code from this project is subject to the terms of the MIT
license. A copy of the MIT license is supplied with the project, or
can be obtained at https://opensource.org/licenses/MIT.
