# metaMAP - quantifying genomes in metagenomes

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow to map and quantify prevalence and abundance of genomes in metagenomes.

Paired metagenomic reads are mapped to a genome database with [BWA](http://bio-bwa.sourceforge.net/) and filtered based on nucleotide identity, mapping score and alignment fraction with [SAMtools](http://www.htslib.org/). Read counts, coverage and depth are calculated per genome for each metagenome and final summary files are generated across samples. Some of these concepts are discussed in great detail by Matt Olm (developer of [dRep](https://drep.readthedocs.io/en/latest/) and [inStrain](https://instrain.readthedocs.io/en/latest/index.html)) [here](https://instrain.readthedocs.io/en/latest/important_concepts.html#an-overview-of-instrain-and-the-data-it-generates).

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

2. Clone repository
```
git clone https://github.com/alexmsalmeida/metamap.git
```

## How to run

1. Edit `config.yml` file to point to the <b>input</b>, <b>output</b> and <b>database</b> directories/files. Input file should contain tab-separated columns with the paths to the forward and reverse reads of each metagenome to map. The <b>database</b> folder should point to a multi-FASTA file containing the reference genomes, indexed with [BWA](http://bio-bwa.sourceforge.net/) (`bwa index`). If mapping to a reference database where each sequence represents its own genome (i.e, when mapping against a viral sequence catalog) change the option `type` in the `config.yml` to "complete".

2. (option 1) Run the pipeline locally (adjust `-j` based on the number of available cores)
```
snakemake --use-conda -k -j 4
```
2. (option 2) Run the pipeline on a cluster (e.g., LSF)
```
snakemake --use-conda -k -j 100 --cluster-config cluster.yml --cluster 'bsub -n {cluster.nCPU} -M {cluster.mem} -o {cluster.output}'
```

## Output

The main output is located in the directory `summary/` which contains four files:

* `bwa_counts-total.csv`: read counts per genome across all samples including multi-mapped reads.
* `bwa_counts-unique.csv`: read counts per genome across all samples excluding multi-mapped reads.
* `bwa_cov-est.csv`: breadth of coverage per genome across all samples.
* `bwa_cov-exp.csv`: expected breadth of coverage per genome across all samples based on their level of read depth. This calculation is further discussed [here](https://instrain.readthedocs.io/en/latest/important_concepts.html#detecting-organisms-in-metagenomic-data).
