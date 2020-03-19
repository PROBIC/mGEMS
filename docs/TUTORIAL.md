# Genomic epidemiology with mixed samples: a tutorial
This tutorial constains instructions on how to repproduce the results
of the three main synthetic experiments from *Genomic epidemiology
with mixed samples*, Mäklin et al. 2020, in preparation.

## Requirements
### mGEMS pipeline
- [Themisto](https://github.com/jnalanko/Themisto)
- [mSWEEP](https://github.com/probic/mSWEEP)
- [mGEMS](https://github.com/probic/mGEMS)
- [shovill](https://github.com/tseemann/shovill/)

### Phylogenetic analysis
- [snippy](https://github.com/tseemann/snippy/)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)

## Tutorial
### Table of Contents

- [Reference data](#referencedata)
- [Synthetic mixed samples](#mixedsamples)
- [Indexing](#indexing)
- [Pseudoalignment](#pseudoalignment)
- [Abundance estimation](#estimation)
- [Binning](#binning)
- [Assembly](#assembly)
- [SNP calling](#snpcalling)
- [Phylogenetic inference](#phylogenetics)


### <a name="referencedata"></a>Reference data
Alternatively, download the reference data from **figshare (TODO:: ADD)**
- [*Escherichia coli*]()
- [*Enterococcus faecalis*]()
- [*Staphylococcus aureus*]()

### <a name="mixedsamples"></a>Synthetic mixed samples
Download the reads from the ENA using the 'get_read.sh' script in the
scripts folder
```
scripts/get_read.sh ERR434377
scripts/get_read.sh ERR435484
scripts/get_read.sh ERR905810
```
Concatenate the read files together
```
zcat ERR434377_1.fastq.gz >> ecoli-1_1.fastq
zcat ERR434377_2.fastq.gz >> ecoli-1_2.fastq
zcat ERR435484_1.fastq.gz >> ecoli-1_1.fastq
zcat ERR435484_2.fastq.gz >> ecoli-1_2.fastq
zcat ERR905810_1.fastq.gz >> ecoli-1_1.fastq
zcat ERR905810_2.fastq.gz >> ecoli-1_2.fastq
```
Compress the mixed sample
```
gzip ecoli-1_1.fastq
gzip ecoli-1_2.fastq
```
Remove the original files
```
rm ERR434377_1.fastq.gz
rm ERR434377_2.fastq.gz
rm ERR435484_1.fastq.gz
rm ERR435484_2.fastq.gz
rm ERR905810_1.fastq.gz
rm ERR905810_2.fastq.gz
```

### <a name="indexing"></a>Indexing
Create a *31*-mer pseudoalignment index with Themisto using four
threads and maximum 8192 megabytes of RAM.
```
mkdir mGEMS-ecoli-reference
mkdir mGEMS-ecoli-reference/tmp
build_index --k 31 --input-file mGEMS-ecoli-reference-v1.0.0.fasta --auto-colors --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp --mem-megas 8192 --n-threads 4
```

### <a name="pseudoalignment"></a>Pseudoalignment
Align the mixed sample files against the index using four threads
```
pseudoalign --query-file ecoli-1_1.fastq.gz --outfile ecoli-1_1.aln --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp -n-threads 4
pseudoalign --query-file ecoli-1_2.fastq.gz --outfile ecoli-1_2.aln --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp --n-threads 4 
```

### <a name="estimation"></a>Abundance estimation
Estimate the relative abundances with mSWEEP and write the results and posterior
probabilities using four threads
```
mSWEEP --themisto-1 ecoli-1_1.aln --themisto-2 ecoli-1_2.aln --themisto-index mGEMS-ecoli-reference -i mGEMS-ecoli-reference-v1.0.0.grouping -o ecoli-1 --write-probs --gzip-probs -t 4
```

### <a name="binning"></a>Binning
Bin the reads with mGEMS and write the binned samples to the
'ecoli-1' folder.
```
mkdir ecoli-1
mGEMS -r ecoli-1_1.fastq.gz,ecoli-1_2.fastq.gz --themisto-alns ecoli-1_1.aln,ecoli-1_2.aln -o ecoli-1 --probs ecoli-1_probs.csv.gz -a ecoli-1_abundances.txt --index mGEMS-ecoli-reference
```
Note this will create bins for **all** reference lineages. If know
which lineages the samples originate from (in our case Escherichia
coli ST131 A, B, and C2), use the '--groups' option to only create
those bins
```
mGEMS --groups Escherichia_coli_ST131-A,Escherichia_coli_ST131-B,Escherichia_coli_ST131-C2 -r ecoli-1_1.fastq.gz,ecoli-1_2.fastq.gz --themisto-alns ecoli-1_1.aln,ecoli-1_2.aln -o ecoli-1 --probs ecoli-1_probs.csv.gz -a ecoli-1_abundances.txt --index mGEMS-ecoli-reference
```

### <a name="assembly"></a>Assembly
Assemble the sequences with shovill using 4 threads and maximum of
8192 megabytes of RAM
```
shovill --outdir ecoli-1/Escherichia_coli_ST131-A --R1 ecoli-1/Escherichia_coli_ST131-A_1.fastq.gz --R2 ecoli-1/Escherichia_coli_ST131-A_2.fastq.gz --cpus 4 --ram 8
shovill --outdir ecoli-1/Escherichia_coli_ST131-B --R1 ecoli-1/Escherichia_coli_ST131-B_1.fastq.gz --R2 ecoli-1/Escherichia_coli_ST131-B_2.fastq.gz --cpus 4 --ram 8
shovill --outdir ecoli-1/Escherichia_coli_ST131-C2 --R1 ecoli-1/Escherichia_coli_ST131-C2_1.fastq.gz --R2 ecoli-1/Escherichia_coli_ST131-C2_2.fastq.gz --cpus 4 --ram 8
```

### <a name="snpcalling"></a>SNP calling
Download the reference sequence 'NCTC13441' from the ENA
```
wget -O NCTC13441.fasta.gz http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/uf/UFZF01.fasta.gz
```
Call SNPs in the genome with snippy
```
mkdir snippy-tmp
snippy --outdir ecoli-1-Escherichia_coli_ST131-A --ctgs ecoli-1/Escherichia_coli_ST131-A/contigs.fa --ref NCTC13441.fasta.gz --cpus 4 --ram 8 --tmpdir snippy-tmp
snippy --outdir ecoli-1-Escherichia_coli_ST131-B --ctgs ecoli-1/Escherichia_coli_ST131-B/contigs.fa --ref NCTC13441.fasta.gz --cpus 4 --ram 8 --tmpdir snippy-tmp
snippy --outdir ecoli-1-Escherichia_coli_ST131-C2 --ctgs ecoli-1/Escherichia_coli_ST131-C@/contigs.fa --ref NCTC13441.fasta.gz --cpus 4 --ram 8 --tmpdir snippy-tmp
```
Build the core SNP alignment with snippy
```
snippy-core --ref ecoli-1-Escherichia_coli_ST131-A/ref.fa ecoli-1-Escherichia_coli_ST131-A ecoli-1-Escherichia_coli_ST131-B ecoli-1-Escherichia_coli_ST131-C2
```

### <a name="phylogenetics"></a>Phylogenetic inference
Use RAxML-NG to infer a maximum likelihood phylogeny from 100 starting
parsimony and random trees under the GTR+G4 model
```
raxml-ng --search --msa core.full.aln --prefix CT --threads 4 --tree rand{100},pars{100} --model GTR+G4
```
Calculate bootstrap support values with 1000 replicates
```
raxml-ng --bootstrap --msa core.full.aln.raxml.rba --bs-trees 1000 --prefix CB --threads 4
```
Perform a bootstrap convergence check
```
raxml-ng --bsconverge --bs-trees CB.raxml.bootstraps --prefix CS --seed 2 --threads 4
```
Add the bootstrap support values to the maximum likelihood tree with
the best likelihood
```
raxml-ng --support --tree CT.raxml.bestTree --bs-trees CB.raxml.bootstraps --prefix CS --threads 4
```
The best tree with the bootstrap values will be written in the
'CS.support' file.
