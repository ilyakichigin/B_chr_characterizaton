# Introduction

DOPseq_analyzer is a set of tools for processing high-throughput sequencing data generated from isolated (flow sorted or microdissected) chromosomes.

Currently, three pipelines are implemented: 

1. Chromosomal region prediction pipeline (`dopseq_pipeline`) includes read trimming, alignment to reference genome, contamination filtering, region calling and statistics calculation.
2. Variant calling and annotation pipeline for validated chromosome-specific regions (`variation_pipeline`).
3. Analysis of anole microchromosomes (`anolis/pipeline/anolis_dopseq_pipe.py`) includes steps similar to dopseq_pipeline. It is optimized for reference genomes consisting of scaffolds and has a possibility to handle WGA libraries. This pipeline is not included in the pipeline installation and can be called only directly. Maintained by ilyakichigin.

# Installation

Dependencies for `dopseq_pipeline`:

1. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested on v.2.1.0, 2.2.4) or [bwa](https://sourceforge.net/projects/bio-bwa/files/) (tested on v.0.7.12)

2. [bedtools](http://bedtools.readthedocs.io/en/latest/) (tested on v.2.17.0 and v.2.24.0)

3. [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) R (Bioconductor) package 

Dependencies for `variation_pipeline`:

1. [GATK 4](https://software.broadinstitute.org/gatk/download/beta) (tested on beta.1)
 
2. [snpEff](http://snpeff.sourceforge.net/) (tested on v.3.3.0 and v.4.3t)

Installation:

```
git clone https://github.com/ilyakichigin/DOPseq_analyzer.git
cd DOPseq_analyzer
pip install --user .
```

# Usage 

Pipeline for chromosome region identification can be called with the command:
```
dopseq_pipeline [-c|-d|-s] dopseq_makefile.yaml
```
Makefile example can be found at `examples/dopseq_makefile.yaml` in this repository or copied to your working directory by running `dopseq_pipeline -c my_makefile.yaml`. Use this file to specify your input data and parameters. 

Use `dopseq_pipeline -s dopseq_makefile.yaml > genome.sizes` to generate tab-separated file listing chromosomes and their sizes for the reference genome specified in the makefile. 

Dry run `-d` option provides command listing and checks if all the input files are present. After that, pipeline can be run with `dopseq_pipeline my_makefile.yaml`

Briefly, this pipeline trims and aligns reads to the reference genome (and optonally contaminant genome, human being most obvious for mammalian chromosome samples), filters the alignment, and classifies selected chromosomes of the reference genome based on mean distances between mapped read positions. Regions with lower means can be further interpreted as present on the isolated chromosomes. Note that these regions cannot be used 'as is' and require manual inspection and correction. Pipeline steps and output files are described in the example makefile. 


Pipeline for variant calling and annotation can be called similarly:
```
variation_pipeline [-d|-c] variation_makefile.yaml
```
Currently, this pipeline runs GATK HaplotypeCaller and snpEff annotation with default settings for a given set of genome regions (presumably present on the chromosome of interest). It also creates a file with variation summary, including per region and per gene statistics. 