# Introduction

DOPseq_analyzer is a a set of tools for processing of sequencing data of isolated (flow sorted or microdissected) chromosomes amplified with Degenerate-Oligonucleotide-Primed PCR (DOP-PCR) or WGA. It divides the reference genome into regions based on density of positions covered by reads. Regions with higher position density can be further interpreted as target (specific to isolated chromosomes). Note that these regions cannot be used 'as is' and require manual inspection and correction. 

Currently, three pipelines are implemented: 

1. Analysis of B chromosomes (pipeline/b_dopseq_pipe.py) includes read trimming, alignment to reference genome, contamination filtering, region calling and statistics calculation. Details on processing steps and additional routines for repetitive DNA characterization with RepeatExplorer and variant calling are described at the end of the configuration file (pipeline/b_dopseq_pipe.config).
2. Analysis of anole microchromosomes (anolis/pipeline/anolis_dopseq_pipe.py) includes similar steps. It is optimized for reference genomes consisting of scaffolds and has a possibility to handle WGA libraries. Maintained by ilyakichigin.
3. Variant calling for validated chromosome-specific regions with GATK HaplotypeCaller (pipeline/vca_reg.py).

# Installation and usage

Dependancies:

1. cutadapt (tested on v.1.8.3)

2. bowtie2 (tested on v.2.1.0, 2.2.4)

3. Pysam python package (tested on v.0.8.1)

4. bedtools (tested on v.2.17.0 and v.2.24.0)

5. DNAcopy R (Bioconductor) package 

6. (Optional - for subsequent variant calling and annotation) GATK (tested on v.3.3.0), picard-tools (tested on v.1.125), snpEff (tested on v.3.3.0)

Scripts located into exec folder can be run independenly or as a part of the pipeline. In the latter case, package structure must be preserved in order for pipeline to find exec folder. Some of the scripts are not included in the pipelines, as these should use corrected regions as input or prepare the data for analysis with outside tools. For more thorough description, see pipeline config files.
