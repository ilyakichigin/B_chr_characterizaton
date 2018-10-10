# Introduction

DOPseq_analyzer is a set of tools for processing high-throughput sequencing data generated from isolated (flow sorted or microdissected) chromosomes.

Current version implements chromosomal region prediction pipeline (`dopseq`) including 
- read trimming (`cutadapt`) and qc (`fastqc`)
- alignment to reference genome (`bwa mem`)
- PCR duplicate and quality filtering (`Picard`, `samtools`), and finally
- region calling with a custom script based on `DNAcopy` Bioconductor package.

This software relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow control and on [conda](https://conda.io/docs/) for dependencies management. The pipeline implementation is based on [Snakemake workflow for dna-seq](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling).

# Quick start

First, install snakemake using your preferred version of [conda](https://conda.io/docs/user-guide/install/index.html)
```
conda install snakemake
```

Clone dopseq, note that the results are written into the installation folder by default
```
git clone https://github.com/ilyakichigin/DOPseq_analyzer.git
cd DOPseq_analyzer
```

Next add your sample data to `samples.tsv` and your analysis parameters to `config.yaml`. After that, you can test the pipeline with a dry run:
```
snakemake --use-conda -n
```

And then run the analysis with
```
snakemake --use-conda
```

The `--use-conda` parameter ensures that for each step all software internally used by `dopseq` will be downloaded and installed inside an isolated environment. Thus, the first run will take longer - in order to get all dependencies in place.

# Steps and outputs description

All output files are stored within the `results` folder, with subfolders numbered according to the analysis order and named by the output type. File names have prefixes corresponding to sample IDs. 

- `0_fastqc_init` - FastQC analysis of the input reads;
- `1_trimmed` - read trimming and filtering with cutadapt;
- `2_fastqc_trim` - FastQC analysis of the trimmed reads, helps to identify the remaining problems with reads;
- `3_mapped` - read mapping to reference genome with bwa mem;
- `4_dedup` - removal of PCR duplicates with Picard MarkDuplicates (BAM removed, stats remain);
- `5_filtered` - final alignment filtering by mapping quality and alignment length;
- `6_merged` - merging all BAMs per sample with samtools merge;
- `7_positions` - overlapping reads into read positions with pybedtools;
- `8_regions` - genome segmentation based on distances between read positions with DNAcopy.

# Output interpretation 

The key outputs of the pipeline are per-sample genome segmentations located in `results/8_regions/{sample}.tsv` files. The segmentation is based on the end-to-start distances between genome positions covered by reads, or PD. The mean values are represented in `pd_mean` column. Regions with lower mean PD are more likely to be present on sampled chromosome (i.e., target regions), while regions with higher mean PD tend to represent a noise resulting from, e.g. contamination with whole-genome or external DNA. Thus, problem of target region identification can transformed into problem of drawing the borderline between low-PD target regions and high-PD contamination. Currently, this step is not automated as there are multiple confounding factors:

- spurious high-PD clusters with low number of markers should be filtered;
- some target regions are over-segmented due to mapping efficiency varying across the chromosome;
- for fragmented reference genome assemblies there is no obvious distinction between high-PD and low-PD regions;
- same problem arises with increase of evolutionary distance between sampled and reference species;
- additional filters based on mean size or mean coverage for mapped positions can be useful in some cases.

# Parameter setting

## samples.tsv

Tab separated file with sample data

`sample` - sample name used as prefix for the output files.

`unit` - multiple inputs (lanes, libraries, biological replicates etc) can be specified per sample. These are trimmed, aligned and filtered separately, so file prefixes for steps 0-5 include both sample and unit names separated by '-'.

`platform` - used to `@RG PL` field of the per-unit BAMs, merged BAM contains all `@RG` lines. 

`adapters` - sequencing or library adapters to be removed from reads:
- `dop` DOP-PCR MW6 primer, reads without primer match at 5\` end are discarded;
- `wga` WGA1 GP primer, reads without primer match at 5\` end are discarded;
- `illumina` only remove Illumina standard adapter from 3\` end;
- `none` do not trim adapters, just filter reads using parameters from configuration file
For paired-end reads trimming with `dop` and `wga`, only pairs with primer matches in both reads are retained. To increase the amount of retained reads, you can specify forward and reverse reads as separate lanes of one sample.

`fq1` - forward or single-end reads fastq file (can be gzipped).

`fq2` - reverse reads fastq file. Keep blank for single-end reads.

## config.yaml

`samples` - path to tab-separated file with sample data

`genome` - path to unpacked reference genome in fasta format. For non-model species selection of the reference balances between evolutionary proximity to the sample species and assembly quality. You may want to experiment with various references in order to obtain better quality results.

`rmdup` - whether to perform PCR duplicate removal, boolean. In theory, should provide more sensible coverage results. However, amplicon-like recovery of same genomic regions in different biological samples suggest the opposite.

`params` section provides software-related parameters. Note that some steps currently use only default settings (fastqc, bwa mem).
- `threads` - number of parallel threads for `bwa mem` and `samtools merge` processes. Total number of cores used by the pipeline is controlled by `snakemake -j` parameter.
- `cutadapt` - cutadapt general filtering options for paired-end (PE) and single-end (SE) reads. Default - trim terminal Ns and remove reads shorter than 20 bp.
- `picard` - remove duplicates instead of marking them
- `filter` - final unit BAM cleanup prior to merging includes filters of minimum mapping quality (use higher value for more stringent removal of repetitive mappings) and mapping length (higher value help to avoid spurious mappings at longer evolutionary distances).

## Other files

This section can be useful if you want to change parameters not listed above, as well as to add or remove steps.

`Snakemake` file sets the desired output files and links to the other smk files: 
- `rules/dopseq.smk` specifies all steps for outputs creation, configurable parameters are automatically picked up from `config.yaml`, for each step conda environment with all dependencies is created either by wrappers (which link to [wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/)), or by conda directive (environments are located in `env` subfolder). 
- `rules/common.smk` contains rules for filename and parameter setting based on sample data, which is located in `samples.tsv`.

Scripts included in the package: `script/regions.py` converts filtered BAM to BED with positions (`results/7_positions`) and performs genome segmentation (`results/8_regions`).


For further reading, please address [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)







```