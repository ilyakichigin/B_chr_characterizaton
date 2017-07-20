#!/usr/bin/env python

# Call and annotate variants for specified genomic regions
# Inputs, options and outputs are listed in the configuration .yaml file
# Requires: GATK (tested on v.3.3.0), picard-tools (tested on v.1.125), samtools
# Usage: vca_reg.py vca_reg.yaml

import os
import sys
import yaml
import shlex
import argparse
import subprocess

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Pipeline for variant calling, filtering, annotation, and alternate reference generation f
                    specified genome regions. 
                    """
                    )
    parser.add_argument("config_file", help="input configuration file")
    parser.add_argument("-d", "--dry_run", action="store_true", default=False,
                        help="Check all dependencies and print out all commands")
    
    return parser.parse_args()

def run_command(command, return_stdout = False, dry_run = True, verbose = True):
    if verbose:
        print command
    command = shlex.split(command)
    if not dry_run:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
        else:
            with open('vca.log','a') as f:
                f.write(err)
        if return_stdout:
            return out

def add_read_groups(bam, gatk, samtools, dry_run = True):
    header = run_command(samtools + ' view -H ' + bam, return_stdout = True, dry_run = False, verbose = False)
    if '@RG' not in header:
        bam_prefix = bam.split('.')[0]
        temp_bam = bam + '.rg.tmp'
        run_command('%s AddOrReplaceReadGroups --RGLB 1 --RGPL Illumina --RGPU 1 --RGSM %s -I %s -O %s'  
            % (gatk, bam_prefix, bam, temp_bam), dry_run = dry_run)
        print 'mv %s %s' % (temp_bam, bam)
        if not dry_run:
            os.rename(temp_bam, bam)

def generate_names(aln, ref):

    names = dict()
    assert aln.endswith('.filter.bam')
    names['pos'] = aln[:-11] + '.pos.bed'
    names['fadict'] = '.'.join(ref.split('.')[:-1]) + '.dict'
    names['md_bam'] = aln[:-4] + '.rmdup.bam'
    names['md_log'] = aln[:-4] + '.rmdup.log'
    names['vcf'] = aln[:-11] + '.reg.vcf'
    names['ann_vcf'] = aln[:-11] + '.reg.ann.vcf'
    names['ann_html'] = aln[:-11] + '.reg.ann.html'
    return names

# parse config

args = parse_command_line_arguments()

with open(args.config_file) as conf_file:
    params = yaml.safe_load(conf_file)
aln = params['bam_file']
ref = params['ref_genome']
reg = params['reg_bed']

fnames = generate_names(aln, ref)

### prepare
# create sequence dictionary for genome
if not os.path.isfile(fnames['fadict']):
    run_command('%s CreateSequenceDictionary -R %s -O %s' 
        % (params['gatk_path'], ref, fnames['fadict']), 
        dry_run = args.dry_run)

# index genome
if not os.path.isfile(ref + '.fai'):
        run_command(params['samtools_path'] + ' faidx ' + ref, dry_run = args.dry_run)

# mark duplicates

if not os.path.isfile(fnames['md_bam']):
    add_read_groups(aln, params['gatk_path'], params['samtools_path'], dry_run = args.dry_run)
    run_command('%s MarkDuplicates -I %s -O %s -M %s'
     % (params['gatk_path'], aln, fnames['md_bam'], fnames['md_log']), dry_run = args.dry_run)

# index dedulpicated BAM file
if not os.path.isfile(fnames['md_bam'] + '.bai'):
    run_command(params['samtools_path'] + ' index ' + fnames['md_bam'], dry_run = args.dry_run)

###analyze
# infer genotypes for verified regions only
if not os.path.isfile(fnames['vcf']):
    run_command('%s HaplotypeCaller -I %s -O %s -R %s -L %s' 
        % (params['gatk_path'], fnames['md_bam'], fnames['vcf'], ref, reg), dry_run = args.dry_run)

# annotate variants
if not os.path.isfile(fnames['ann_vcf']):
    ann_data = run_command('java -Xmx%s -jar %s -s %s %s %s ' 
        % (params['max_mem'], params['snpEff_path'], fnames['ann_html'], params['snpEff_genome'], fnames['vcf']), 
        return_stdout = True, dry_run = args.dry_run)
    if not args.dry_run:
        with open(fnames['ann_vcf'], 'w') as out:
            print 'Writing to ' + fnames['ann_vcf']
            out.write(ann_data)