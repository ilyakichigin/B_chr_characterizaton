# Introduction

DOPseq_analyzer is a a set of tools for processing of Illumina sequencing data obtained from Degenerate Oligonucleotide Primed PCR libraries of isolated (flow sorted or microdissected) chromosomes. The result is the division of the genome used as reference into regions, which can be further interpreted as target (specific to isolated chromosomes) or contamination (non-specific). Note that these regions cannot be used 'as is' and require manual inspection and correction.

Currently, two pipelines are implemented: 

1. Analysis of B chromosomes (pipeline/b_dopseq_pipe.py) includes read trimming, alignment to reference genome, contamination filtering, region calling and statistics calculation. It is suitable for reference genomes assembled up to the chromosomes.
2. Analysis of anole microchromosomes (anolis/pipeline/anolis_dopseq_pipe.py) includes similar steps. It is optimized for reference genomes consisting of scaffolds and has a possibility to handle WGA libraries.

# Installation and usage

Dependancies:

1. cutadapt (tested on v.1.6)

2. bowtie2 (tested on v.2.1.0, 2.2.4)

3. Pysam python package (tested on v.0.8.1)

4. bedtools (tested on v.2.17.0)

5. DNAcopy R package 

6. (Optional - for subsequent variant calling) GATK (tested on v.3.3.0) and picard-tools (tested on v.1.125)

Scripts located into exec folder can be run independenly or as a part of the pipeline. Some of the scripts are not included in the pipelines, as these should use corrected regions as input. For more thorough description, see pipeline config files.
