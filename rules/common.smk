import pandas as pd
from snakemake.utils import validate

# report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
# reactivate later with meaningful parameters
# validate(config, schema="../schemas/config.schema.yaml")

# samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["samples"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
units['prefix'] = units['sample'] + '-' + units['unit']
validate(units, schema="../schemas/samples.schema.yaml")

# contigs in reference genome
# contigs = pd.read_table(config["genome"] + ".fai",
#                         header=None, usecols=[0], squeeze=True, dtype=str)


##### Wildcard constraints #####
wildcard_constraints:
    # vartype="snvs|indels",
    sample="|".join(units["sample"].unique()),
    unit="|".join(units["unit"]),
    # contig="|".join(contigs)


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}-{unit}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, unit=wildcards.unit,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("results/1_trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "results/1_trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


# def get_dedup_bams(wildcards):
#     """Get all aligned reads of given sample."""
#     return expand("results/4_dedup/{sample}-{unit}.bam",
#                   sample=wildcards.sample,
#                   unit=units.loc[wildcards.sample].unit)

def get_filtered_bams(wildcards):
    """Get all per-unit alignments of given sample"""
    return expand("results/5_filtered/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)

def get_position_beds(wildcards):
    """Get all position BED files"""
    return expand("results/7_positions/{sample}.bed",
                 sample=units["sample"].unique())

def ampl_to_cutadapt_pe(wildcards):

    platform=units.loc[(wildcards.sample, wildcards.unit), "platform"]
    ampl=units.loc[(wildcards.sample, wildcards.unit), "adapters"]
    # anchored linked adapters
    if ampl == 'dop':
        return '-a CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG -A CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG --discard-untrimmed'
    elif ampl == 'wga':
        return '-a TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -A TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -e 0.2 --discard-untrimmed'
    # relaxed versions - do not discard untrimmed
    # elif ampl == 'dop_relaxed':
    #     return '-a CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG -A CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG'
    # elif ampl == 'wga_relaxed':
    #     return '-a TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -A TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -e 0.2'
    elif ampl == 'illumina':
        return ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    return ''

def ampl_to_cutadapt_se(wildcards):

    platform=units.loc[(wildcards.sample, wildcards.unit), "platform"]
    ampl=units.loc[(wildcards.sample, wildcards.unit), "adapters"]
    # anchored linked adapters
    if ampl == 'dop':
        return '-a CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG --discard-untrimmed'
    elif ampl == 'wga':
        return '-a TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -e 0.2 --discard-untrimmed'
    # relaxed versions - do not discard untrimmed
    # elif ampl == 'dop_relaxed':
    #     return '-g CCGACTCGAGNNNNNNATGTGG...CCACATNNNNNNCTCGAGTCGG'
    # elif ampl == 'wga_relaxed':
    #     return '-g TTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACAA -e 0.2'
    elif ampl == 'illumina':
        return ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    return ''
