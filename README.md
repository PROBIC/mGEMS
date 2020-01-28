# msweep-assembly

mSWEEP genome assembly plugin code.

# Installation
## Compiling from source
### Requirements
- C++11 compliant compiler.
- cmake

### Compilation
Clone the repository
```
git clone https://github.com/PROBIC/msweep-assembly.git
```
enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
This will compile the read_alignment, assign_reads, and build_sample executables in the build/bin/ directory.


# Usage
Align paired-end reads 'reads_1.fastq.gz' and 'reads_2.fastq.gz' with [Themisto]()
```
pseudoalign --index-dir themisto_index --query-file reads_1.fastq.gz --outfile pseudoalignments_1.txt --rc --temp-dir tmp --n-threads 16 --mem-megas 8192
pseudoalign --index-dir themisto_index --query-file reads_2.fastq.gz --outfile pseudoalignments_2.txt --rc --temp-dir tmp --n-threads 16 --mem-megas 8192
```

Convert the pseudoalignment to [kallisto]() format using [telescope]()
```
ntargets=$(sort themisto_index/coloring-names.txt | uniq | wc -l)
telescope --n-refs $ntargets -r pseudoalignments_1.txt,pseudoalignments_2.txt -o outfolder --mode intersection
```

Create a fake kallisto-style run_info.json file
```
Themisto_run_info.sh $(wc -l outfolder_1.txt) $ntargets > outfolder/run_info.json
```

Determine read assignments to equivalence classes from the kallisto
format files
```
read_alignment -e outfolder/outfolder.ec -s outfolder/read-to-ref.txt -o outfolder --write-ecs --themisto --n-refs $ntargets --gzip-output
```

Estimate the relative abundances with mSWEEP
```
mSWEEP -f outfolder -i reference_grouping.txt -o msweep-out --write-probs --gzip-probs
```

Extract the names of the 3 most abundant reference groups
```
grep -v "^[#]" msweep-out_abundances.txt | sort -rgk2 | cut -f1 | head -n3 > most_abundant_groups.txt
```

Assign reads to the 3 most abundant reference groups based on the estimated probabilities
```
assign_reads -f outfolder/ec_to_read.csv.gz -p msweep-out_probs.csv.gz -a msweep-out_abundances.txt -o outfolder/ --groups most_abundant_groups.txt --gzip-output
```

Construct the binned samples from the original files

```
while read -r sample; do
	build_sample -a outfolder/$sample\"\"_reads.txt.gz -o outfolder/$sample -1 reads_1.fastq.gz -2 reads_2.fastq.gz --gzip-output
done < most_abundant_groups.txt
```
