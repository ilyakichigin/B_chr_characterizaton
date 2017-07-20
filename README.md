# Introduction

DOPseq_analyzer is a a set of tools for processing of sequencing data of isolated (flow sorted or microdissected) chromosomes amplified with Degenerate Oligonucleotide Primed PCR (DOP-PCR) or WGA.

Currently, three pipelines are implemented: 

1. Analysis of B chromosomes (pipeline/b_dopseq_pipe.py) includes read trimming, alignment to reference genome, contamination filtering, region calling and statistics calculation.
2. Analysis of anole microchromosomes (anolis/pipeline/anolis_dopseq_pipe.py) includes similar steps. It is optimized for reference genomes consisting of scaffolds and has a possibility to handle WGA libraries. Maintained by ilyakichigin.
3. Variant calling for validated chromosome-specific regions with GATK HaplotypeCaller (pipeline/vca_reg.py).

# Installation

Dependencies:

1. [cutadapt](http://cutadapt.readthedocs.io/en/stable/) (tested on v.1.8.3)

2. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (tested on v.2.1.0, 2.2.4)

3. [pysam](http://pysam.readthedocs.io/en/latest/api.html) python package (tested on v.0.8.1)

4. [bedtools](http://bedtools.readthedocs.io/en/latest/) (tested on v.2.17.0 and v.2.24.0)

5. [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) R (Bioconductor) package 

6. (Optional - for subsequent variant calling and annotation) [GATK 4](https://software.broadinstitute.org/gatk/download/beta) (tested on beta.1), [snpEff](http://snpeff.sourceforge.net/) (tested on v.3.3.0)


After dependencies are installed, individual scrips and pipelines can be used. Note that internal folder structure must be preserved so that b_dopseq_pipe.py could find executable scripts. 

# Usage 

Pipeline for chromosome region identification can be called with the command:
```
/some/folder/DOPseq_analyzer/pipeline/b_dopseq_pipe.py [-d] b_dopseq_pipe.yaml
```
Details on inputs, outputs and assignable parameters can be found in example `b_dopseq_pipe.yaml` config file.   Enabling dry run `-d` option provides command listing and checks if all the input files are present. 

Briefly, this pipeline trims and aligns reads to the reference genome (and optonally potential contaminant genome), filters the alignment, and classifies selected chromosomes of the reference genome based on mean distances between mapped read positions. Regions with lower means can be further interpreted as specific to the isolated chromosomes. Note that these regions cannot be used 'as is' and require manual inspection and correction.


Pipeline for variant calling and annotation can be called similarly:
```
/some/folder/DOPseq_analyzer/pipeline/vca_reg.py [-d] vca_reg.yaml
```
Currently, this pipeline simply runs GATK HaplotypeCaller and snpEff annotation with default settings for a given set of regions (presumably present on the chromosome of interest). 
