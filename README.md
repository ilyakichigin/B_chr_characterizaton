# Introduction

DOPseq_analyzer is a set of tools for processing high-throughput sequencing data generated from isolated (flow sorted or microdissected) chromosomes.

Current version implements chromosomal region prediction pipeline: 
- read trimming and qc,
- alignment to reference genome,
- PCR duplicate and quality filtering, and finally
- genome segmentation.

This software relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow control and on [conda](https://conda.io/docs/) for dependencies management. The pipeline implementation is based on [Snakemake workflow for dna-seq](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling).

# Quick start

First, install snakemake using your preferred version of [conda](https://conda.io/docs/user-guide/install/index.html)
```
conda install -c bioconda -c conda-forge snakemake
```

Clone dopseq, checkout the snakemake version branch.
```
git clone https://github.com/ilyakichigin/DOPseq_analyzer.git
cd DOPseq_analyzer
git checkout smk_dev
```
Note that the results are written into the installation folder by default.

Next, add your sample data to `samples.tsv` and your analysis parameters to `config.yaml`. See ['Parameter setting'](#Parameter_setting) section for details.

After that, you can test the pipeline with a dry run:
```
snakemake --use-conda -n
```

And then run the analysis with
```
snakemake --use-conda
```

`--use-conda` parameter ensures that for each step nececcary packages will be downloaded and installed inside an isolated environment. Thus, the first run will take longer - in order to get all dependencies in place.

# Steps and outputs description

Output files are stored within the `results` folder, with subfolders numbered according to the analysis sequence and named by the output type. File names have prefixes corresponding to sample IDs. 

- `0_fastqc_init` - FastQC on the input reads;
- `1_trimmed` - read trimming and filtering with cutadapt;
- `2_fastqc_trim` - FastQC on the trimmed reads, helps to identify the remaining problems with reads;
- `3_mapped` - read mapping to reference genome with bwa mem;
- `4_dedup` - removal of PCR duplicates with Picard MarkDuplicates (BAM removed, stats remain);
- `5_filtered` - alignment filtering by mapping quality and aligned fragment length with samtools and awk;
- `6_merged` - merging all BAMs per sample with samtools;
- `7_positions` - merging overlapping reads into read positions with pybedtools;
- `8_regions` - genome segmentation based on distances between read positions with DNAcopy;
- `stats.xlsx` - statistics for sequencing, mapping, filtering, and positions.

# Output interpretation 

The key outputs of the pipeline are per-sample genome segmentations located in `results/8_regions/{sample}.tsv` files. The segmentation is based on the end-to-start distances between genome positions covered by reads, or PD. The mean values of PD per segment in `pd_mean` column. Regions with lower mean PD are more likely to be present on sampled chromosome and will be subsequently called _target regions_. Regions with higher mean PD tend to represent a noise that can result from erroneous mappings or contamination with whole-genome or external DNA. Thus, problem of target region identification can transformed into problem of drawing the borderline between low 'pd_mean' target regions and high 'pd_mean' contamination. Currently, this step is not automated as there are multiple confounding factors:

- spurious high-PD segments with low number of markers should be filtered, but threshold varies;
- additional filters based on mean position size or mean position coverage for mapped positions can be useful in some cases;
- segmentation for small contigs/scaffolds is unreliable, threshold of 50 kbp was mostly used, but it cannot be recommended for all cases;
- for highly fragmented reference genome assemblies there can be no obvious distinction between high-PD and low-PD regions;
- same problem arises with increase of evolutionary distance between sampled and reference species;
- some target regions are over-segmented due to mapping efficiency varying across the chromosome and should be joined;
- region margins sometimes should be manually corrected.

General recommendation is to sort region prediction files by 'pd_mean', set up filters that fit your data (min 'num_mark', min 'chrom_len', optionally - min 'pos_len_mean' and min 'prop_reg_covered'), and extract top most significant target regions by looking for significant gap in 'pd_mean' (e.g., several target regions have 1-5 kbp 'pd_mean', then first non-target regions has 20 kbp 'pd_mean', then values increase gradually). For manual correction of region margins, visualization of BAMs (step 6) and BEDs (step 7) in genome browser can be useful. 

# Parameter setting

## samples.tsv

Tab separated file with sample data:

`sample` - sample name used as prefix for the output files.

`unit` - multiple input units (lanes, libraries, biological replicates etc) can be specified for each sample. Units are trimmed, aligned and filtered separately, so file prefixes for steps 0-5 include both sample and unit names separated by '-'. Repeated setting of identical sample-unit combinations results in an obscure "ValueError: The truth value of a Series is ambiguous".

`platform` - used to `@RG PL` field of the per-unit BAMs, merged BAM contains all `@RG` lines. 

`adapters` - sequencing or library adapters to be removed from reads:
- `dop` DOP-PCR MW6 primer, reads without primer match at 5\` end are discarded;
- `wga` WGA1 GP primer, reads without primer match at 5\` end are discarded;
- `illumina` only remove Illumina standard adapter from 3\` end;
- `none` do not trim adapters, just filter reads using parameters from configuration file

For paired-end reads trimming with `dop` and `wga`, only pairs with primer matches in both reads are retained. To increase the amount of retained reads, you can specify forward and reverse reads as separate units of one sample.

`fq1` - forward or single-end reads fastq file (can be plain or gzipped).

`fq2` - reverse reads fastq file. Keep blank for single-end reads.

## config.yaml

`samples` - path to tab-separated file with sample data.

`genome` - path to unpacked reference genome in fasta format. For non-model species selection of the reference balances between evolutionary proximity to the sample species and assembly quality. You may want to experiment with various references in order to obtain better quality results. Only single reference genome per project is currently supported.

`rmdup` - whether to perform PCR duplicate removal, boolean. In theory, should provide more sensible coverage statistics. However, amplicon-like recovery of same genomic regions in different biological samples suggest the opposite.

`params` section provides software-related parameters:
- `threads` - number of parallel threads for `bwa mem` and `samtools merge` processes. Total number of cores used by the pipeline is controlled by snakemake `-j` parameter.
- `cutadapt` - cutadapt general filtering options for paired-end (PE) and single-end (SE) reads. Default - trim terminal Ns and remove reads shorter than 20 bp.
- `picard` - by default set to remove duplicates instead of marking them.
- `filter` - includes filters of minimum mapping quality (use higher value for more stringent removal of repetitive mappings) and mapping length (higher value help to avoid mapping errors at longer evolutionary distances).

## Other files

This section can be useful if you want to change parameters not listed above, as well as to add or remove steps.

`Snakemake` file sets the desired output files and links to the other smk files: 
- `rules/dopseq.smk` specifies all steps for outputs creation, parameters are automatically picked up from `config.yaml`. For each step conda environment with all dependencies is created either by wrapper (which is linked to [wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/)), or by conda section (environments are located in `env` subfolder). 
- `rules/common.smk` contains functions for filename and parameter setting based on sample data, which is located in `samples.tsv`.

Scripts included in the package: 
- `script/regions.py` converts filtered BAM to BED with overlapping reads merged into positions (`results/7_positions`) and performs genome segmentation (`results/8_regions`).
- `script/regions.py` collects statistics per unit and per sample in a single xlsx file.

For further reading, please go to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).







```
