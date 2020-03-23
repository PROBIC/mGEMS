# Genomic epidemiology with mixed samples: a tutorial
This tutorial constains instructions on how to repproduce the results
of the three main synthetic experiments from *Genomic epidemiology
with mixed samples*, Mäklin et al. 2020, in preparation.

The tutorial will focus on reproducing the *Escherichia coli*
experiment but contains instructions on how to adapt the scripts to
the *Enterococcus faecalis* and *Staphylococcus aureus* experiments.

For quick instructions on how to run the pipeline in a general
setting, please refer to the README.md file in the root of this
repository.

## Requirements
### mGEMS pipeline
- [Themisto](https://github.com/jnalanko/Themisto)
- [mSWEEP](https://github.com/probic/mSWEEP)
- [mGEMS](https://github.com/probic/mGEMS)
- [shovill](https://github.com/tseemann/shovill/)

### Phylogenetic analysis
- [snippy](https://github.com/tseemann/snippy/)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)

### Extra (macOS only)
#### GNU coreutils
On a macOS system, you'll also need to install GNU coreutils from
homebrew and alias the macOS zcat command to the GNU zcat command for
the duration of the session
```
brew install coreutils
alias zcat=gzcat
ulimit -n 2048
```
#### Concurrent file connections limit
macOS also limits the number of concurrent file connections, which
will have to be increased to run Themisto and shovill
```
ulimit -n 2048
```

## Tutorial
### Table of Contents

- [Select a species](#selectspecies)
- [Reference data](#referencedata)
- [Synthetic mixed samples](#mixedsamples)
- [Indexing](#indexing)
- [Pseudoalignment](#pseudoalignment)
- [Abundance estimation](#estimation)
- [Binning](#binning)
- [Assembly](#assembly)
- [SNP calling](#snpcalling)
- [Phylogenetic inference](#phylogenetics)


### <a name="selectspecies"></a>Select a species
Download the supplementary table from the mGEMS manucsript which
contains the relevant information
```
https://zenodo.org/record/3724144/files/mGEMS_Supplementary_Table_mixed_samples.tsv
```
Filter the table to contain only the *E. coli* (ecoli) experiments
```
grep "ecoli" mGEMS_Supplementary_Table_mixed_samples.tsv" > mixed_samples.tsv
```
If you want to reproduce the *E. faecalis* experiments, change 'ecoli'
to 'efaec'. For *S. aureus*, change 'ecoli' to 'saur'. Running these
other two experiments may require resources beyond the typical laptop or
desktop computer.

### <a name="referencedata"></a>Reference data
Download and extract the relevant reference data from zenodo
- [*E. coli*](https://zenodo.org/record/3724112)
- [*E. faecalis*](https://zenodo.org/record/3724101)
- [*S. aureus*](https://zenodo.org/record/3724135)

by running
```
wget https://zenodo.org/record/3724112/files/mGEMS-ecoli-reference-v1.0.0.tar.gz
tar -zxvf mGEMS-ecoli-reference-v1.0.0.tar.gz
```
Construction of the reference dataset(s) is describe in more detail in
Mäklin et al. 2020.

### <a name="indexing"></a>Indexing
Create a *31*-mer pseudoalignment index with Themisto using two
threads and maximum 8192 megabytes of RAM.
```
mkdir mGEMS-ecoli-reference
mkdir mGEMS-ecoli-reference/tmp
build_index --k 31 --input-file mGEMS-ecoli-reference-sequences-v1.0.0.fasta.gz --auto-colors --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp --mem-megas 8192 --n-threads 2
```
change 'ecoli' to 'efaec' or 'saur' if you are trying to reproduce the
other experiments.

### <a name="mixedsamples"></a>Synthetic mixed samples
Download the isolate sequencing data and create the synthetic mixed
samples by concatenating the isolate files
```
## Download the sequencing data and create the samples
oldid=""
while read line; do
	id=$(echo $line | cut -f3 -d' ')
	sample=$(echo $line | cut -f1 -d' ')
	scripts/get_forward.sh $sample | gunzip -c >> $id""_1.fastq
	scripts/get_reverse.sh $sample | gunzip -c >> $id""_2.fastq
	if [[ "$id" != "$oldid" ]]; then
		if [ ! -z "$oldid" -a "$oldid" != "" ]; then
			gzip $oldid""_1.fastq
			gzip $oldid""_2.fastq
		fi
	fi
	oldid=$id
done < mixed_samples.tsv
gzip $oldid""_1.fastq
gzip $oldid""_2.fastq
```

### <a name="pseudoalignment"></a>Pseudoalignment
Align the mixed sample files against the index using two threads
```
for f1 in *_1.fastq.gz; do
	f=${f1%_1.fastq.gz}
	f2=$f""_2.fastq.gz
	pseudoalign --query-file $f1 --outfile $f""_1.aln --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp --n-threads 2 --rc
	pseudoalign --query-file $f2 --outfile $f""_2.aln --index-dir mGEMS-ecoli-reference --temp-dir mGEMS-ecoli-reference/tmp --n-threads 2 --rc
	gzip $f""_1.aln
	gzip $f""_2.aln
done
```

### <a name="estimation"></a>Abundance estimation
Estimate the relative abundances with mSWEEP and write the results and posterior
probabilities using two threads
```
for f1 in *_1.fastq.gz; do
	f=${f1%_1.fastq.gz}
	mkdir $f
	mSWEEP --themisto-1 $f""_1.aln.gz --themisto-2 $f""_2.aln.gz --themisto-index mGEMS-ecoli-reference -i mGEMS-ecoli-reference-grouping-v1.0.0.txt -o $f/$f --write-probs --gzip-probs -t 2
done
```

### <a name="binning"></a>Binning
Bin the reads with mGEMS and write the binned samples to the
'ecoli-1' folder.
```
while read line; do
	id=$(echo $line | cut -f3 -d' ')
	cluster=$(echo $line | cut -f2 -d' ')
	mGEMS --groups $cluster -r $id""_1.fastq.gz,$id""_2.fastq.gz --themisto-alns $id""_1.aln.gz,$id""_2.aln.gz -o $id --probs $id/$id""_probs.csv.gz -a $id/$id""_abundances.txt --index mGEMS-ecoli-reference
done < mixed_samples.tsv
```
Note that by default mGEMS creates bins for **all** reference lineages. If know
which lineages the samples originate from (in our case these are
supplied in the mixed_samples.tsv in the second column), the
'--groups' option enables you to only create those bins. Multiple
groups can be binned in the single run by supplying them as a
comma-separated list.

### <a name="assembly"></a>Assembly
Assemble the sequences with shovill using 2 threads and maximum of
8192 megabytes of RAM
```
while read line; do
	id=$(echo $line | cut -f3 -d' ')
	cluster=$(echo $line | cut -f2 -d' ')
	shovill --outdir $id/$cluster --R1 $id/$cluster""_1.fastq.gz --R2 $id/$cluster""_2.fastq.gz --cpus 2 --ram 8
	mv $id/$cluster/contigs.fa $id/
	rm -rf $id/$cluster
	mkdir $id/$cluster
	mv $id/contigs.fa $id/$cluster/
done < mixed_samples.tsv
```

### <a name="snpcalling"></a>SNP calling
Download the reference sequence 'NCTC13441' from the ENA
```
wget -O NCTC13441.fasta.gz http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/uf/UFZF01.fasta.gz
```
Call SNPs in the genome with snippy
```
mkdir snippy-tmp
gunzip NCTC13441.fasta.gz
while read line; do
	id=$(echo $line | cut -f3 -d' ')
	sample=$(echo $line | cut -f1 -d' ')
	cluster=$(echo $line | cut -f2 -d' ')
	snippy --outdir $id/$cluster/$sample --ctgs $id/$cluster/contigs.fa --ref NCTC13441.fasta --cpus 2 --ram 8 --tmpdir snippy-tmp
done < mixed_samples.tsv
gzip NCTC13441.fasta
```
Build the core SNP alignment with snippy
```
snippys=""
while read line; do
	id=$(echo $line | cut -f3 -d' ')
	sample=$(echo $line | cut -f1 -d' ')
	cluster=$(echo $line | cut -f2 -d' ')
	snippys=$snippys""$id/$cluster/$sample" "
	last=$id/$cluster/$sample
done < mixed_samples.tsv
snippy-core --ref $last/ref.fa $snippys
```
the alignment will be stored in the 'core.full.aln' file.

### <a name="phylogenetics"></a>Phylogenetic inference
Use RAxML-NG to infer a maximum likelihood phylogeny from 10 starting
parsimony and random trees under the GTR+G4 model
```
raxml-ng --search --msa core.full.aln --prefix CT --threads 2 --tree rand{10},pars{10} --model GTR+G4
```
Calculate bootstrap support values with 100 replicates
```
raxml-ng --bootstrap --msa core.full.aln --bs-trees 100 --prefix CB --threads 2 --model GTR+G4
```
Perform a bootstrap convergence check
```
raxml-ng --bsconverge --bs-trees CB.raxml.bootstraps --prefix CS --seed 2 --threads 2
```
Add the bootstrap support values to the maximum likelihood tree with
the best likelihood
```
raxml-ng --support --tree CT.raxml.bestTree --bs-trees CB.raxml.bootstraps --prefix CS --threads 2
```
The best tree with the bootstrap values will be written in the
'CS.raxml.support' file.
