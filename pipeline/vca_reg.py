#!/usr/bin/env python

# Call and annotate variants for specified genomic regions
# Inputs, options and outputs are listed in the configuration .yaml file
# Requires: GATK (tested on v.3.3.0), picard-tools (tested on v.1.125), samtools
# Usage: vca_reg.py vca_reg.yaml

import os
import sys
import yaml
import shlex
import subprocess

def run_command(command, return_stdout = False, verbose = True):
    if verbose:
        print ' '.join(command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process.communicate()
    if process.returncode != 0:
        sys.stderr.write(err)
        sys.exit()
    if return_stdout:
        return out

def run_java_program(jarfile, options, max_mem = "1g", return_stdout = False):
    command = ['java', '-Xmx' + max_mem, '-jar', jarfile] + shlex.split(options) # shlex str split for correct quoted string interpretation
    return run_command(command, return_stdout = return_stdout)

def add_read_groups(bam, picard, samtools):
    header = run_command([samtools, 'view', '-H', bam], return_stdout = True, verbose = False)
    if '@RG' not in header:
        bam_prefix = bam.split('.')[0]
        temp_bam = bam + '.rg.tmp'
        run_java_program(picard,
            'AddOrReplaceReadGroups RGLB=1 RGPL=Illumina RGPU=1 RGSM=%s I=%s O=%s' % 
            (bam_prefix, bam, temp_bam))
        print 'mv %s %s' % (temp_bam, bam)
        os.rename(temp_bam, bam)


# parse config
with open(sys.argv[1]) as conf_file:
    params = yaml.safe_load(conf_file)
rg = params['ref_genome']
sample = params['sample']
samtools = params['samtools_path']
picard = params['picard_path']
gatk = params['gatk_path']

# create sequence dictionary for genome
fadict = '.'.join(rg.split('.')[:-1]) + '.dict'
if not os.path.isfile(fadict):
    run_java_program(picard,
        'CreateSequenceDictionary R=%s O=%s' % (rg, fadict))

# index genome
if not os.path.isfile(rg + '.fai'):
        run_command([samtools, 'faidx', rg])


for aln in params['bam_files']:
    
    # add read group data to BAM file
    assert aln.endswith('.bam')
    add_read_groups(aln, picard, samtools)
    
    # mark duplicates in BAM file
    md_bam = aln[:-4] + '.rmdup.bam'
    if not os.path.isfile(md_bam):
        md_log = aln[:-4] + '.rmdup.txt'
        run_java_program(picard,
            'MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s' % (aln, md_bam, md_log))
    
    # index dedulpicated BAM file
    if not os.path.isfile(md_bam + '.bai'):
        run_command([samtools, 'index', md_bam])

    # infer genotypes at bp-resolution for target regions only
    gvcf = aln[:-4] + '.g.vcf'
    if not os.path.isfile(gvcf):
        run_java_program(gatk, 
            '-R %s -T HaplotypeCaller -I %s -ERC BP_RESOLUTION -L %s -o %s' % 
                (rg, md_bam, params['reg_bed'], gvcf),
            max_mem = params['gatk_mem'])


# call variants for all samples into one vcf
vcf = sample + '.vcf'
if not os.path.isfile(vcf):
    gvcfs = ['--variant ' + aln[:-4] + '.gvcf' for aln in params['bam_files']]
    run_java_program(gatk, 
        '-R %s -T GenotypeGVCFs %s -o %s' % 
            (rg, ' '.join(gvcfs), vcf),
        max_mem = params['gatk_mem'])

# filter variants
vcf_filtered = sample + '.hq.vcf'
if not os.path.isfile(vcf_filtered):
    run_java_program(gatk, 
        '-R %s -T SelectVariants %s --variant %s -o %s' % 
            (rg, ' '.join(params['gatk_filters']), vcf, vcf_filtered))

# annotate filtered variants
vcf_annotated = sample + '.hq.ann.vcf'
if not os.path.isfile(vcf_annotated):
    ann_data = run_java_program(params['snpEff_path'],
        '-s %s %s %s ' % (sample + '.hq.ann.summary.html', params['snpEff_genome'], vcf_filtered), 
        return_stdout = True)
    with open(vcf_annotated, 'w') as out:
        print 'Writing to ' + vcf_annotated
        out.write(ann_data)
