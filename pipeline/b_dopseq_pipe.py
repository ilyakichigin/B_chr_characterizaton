#!/usr/bin/env python

import os
import sys
import yaml
import pysam
import argparse
import subprocess

exec_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"..","exec"))
sys.path.append(exec_path)

import fastq_clean
import fastq_aln
import contam_filter
import bam_to_bed
import sample_stats

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Pipeline for processing of sequencing data of DOP-PCR libraries from isolated chromosomes.
                    See config for process description, inputs and outputs. 
                    """
                    )
    parser.add_argument("config_file", help="input configuration file")
    parser.add_argument("-d", "--dry_run", action="store_true", default=False,
                        help="Check all dependencies and print out all commands")
    
    return parser.parse_args()

#def run_module_main(module, args):

def generate_filenames(sample, target_path, contam_path):

    names = dict()
    names['f_trim_fq'] = sample + '.ca.R1.fastq.gz'
    names['r_trim_fq'] = sample + '.ca.R2.fastq.gz'
    names['target_path'] = target_path
    names['contam_path'] = contam_path
    names['target_name'] = target_path.split('/')[-1].split('.')[0]
    names['target_index_prefix'] = '/'.join(target_path.split('/')[:-1]) + names['target_name']
    names['contam_name'] = contam_path.split('/')[-1].split('.')[0]
    names['contam_index_prefix'] = '/'.join(target_path.split('/')[:-1]) + names['contam_name']
    names['target_bam'] = '.'.join([sample,names['target_name'],'bam'])
    names['contam_bam'] = '.'.join([sample,names['contam_name'],'bam'])
    names['filtered_bam'] = '.'.join([sample,names['target_name'],'filter','bam'])
    names['pos_bed'] = '.'.join([sample,names['target_name'],'pos','bed'])
    names['reg_tsv'] = '.'.join([sample,names['target_name'],'reg','tsv'])
    names['reg_pdf'] = '.'.join([sample,names['target_name'],'reg','pdf'])
    names['stat_file'] = conf['sample'] + '.stats.txt'

    return names

if __name__ == '__main__':

    # parse arguments and config
    args = parse_command_line_arguments()
    with open(args.config_file) as conf_file:
        conf = yaml.safe_load(conf_file)

    # dry run notification
    if args.dry_run:
        sys.stderr.write('This is a dry run. Only command listing will be produced.\n')


    # ?add executables check here?
    
    # generate names
    if 'contam_genome' not in conf.keys(): # no contam filtering
        conf['contam_genome'] = conf['target_genome']
    fnames = generate_filenames(conf['sample'],conf['target_genome'],conf['contam_genome'])
    if 'fastq_R_file' not in conf.keys(): # single-end input reads
        conf['fastq_R_file'] = None
        fnames['r_trim_fq'] = None

    if not os.path.isfile(fnames['filtered_bam']):
        sys.stderr.write('----fastq_clean.py----\n') # trim adapters
        fc_args = argparse.Namespace(fastq_F_file=conf['fastq_F_file'], fastq_R_file=conf['fastq_R_file'],
                                     sample_name=conf['sample'], path_to_cutadapt=conf['cutadapt_path'],
                                     ampl=conf['ampl'], params=conf['cutadapt_args'],
                                     dry_run=args.dry_run, trim_illumina=True)
        
        fastq_clean.main(fc_args)
        sys.stderr.write('----Complete!----\n')

        sys.stderr.write('----fastq_aln.py----\n') # align reads
        if conf['aln'] == 'bt2':
            aligner_path = conf['bowtie2_path']
            aligner_args = conf['bowtie2_args']
        elif conf['aln'] == 'bwa':
            aligner_path = conf['bwa_path']
            aligner_args = conf['bwa_args']
        elif conf['aln'] == 'bbm':
            aligner = conf['bbmap_path']
            aligner_args = conf['bbmap_args']
        else:
            raise Exception('Invalid aligner %s!' % conf['aln'])
        for g in ('target', 'contam'):
            fa_args = argparse.Namespace(fastq_F_file=fnames['f_trim_fq'], 
                fastq_R_file=fnames['r_trim_fq'], aligner=conf['aln'], 
                reference_genome=fnames[g + '_path'], path_to_aligner=aligner_path, 
                aligner_args=aligner_args, dry_run = args.dry_run)
            fastq_aln.main(fa_args)
        sys.stderr.write('----Complete!----\n')

        sys.stderr.write('----contam_filter.py----\n') # filter contamination
        cf_args = argparse.Namespace(target_bam=fnames['target_bam'], contam_bam=fnames['contam_bam'],
                                     min_quality=conf['min_mapq'], min_length=conf['min_len'], dry_run = args.dry_run)
        contam_filter.main(cf_args)
        sys.stderr.write('----Complete!----\n')
        '''
        os.remove(fnames['target_sam'])
        if conf['contam_genome'] != conf['target_genome']:
            os.remove(fnames['contam_sam'])
        '''
    else:
        sys.stderr.write('%s filtered alignment to target genome exists. Skipping steps for its generation.\n' % fnames['filtered_bam'])
    sys.stderr.write('----bam_to_bed.py----\n') # create bed file with positions covered by reads
    btb_args = argparse.Namespace(bam_file=fnames['filtered_bam'], bedtools_path='bedtools', 
        dry_run = args.dry_run)
    bam_to_bed.main(btb_args)
    sys.stderr.write('----Complete!----\n')

    # predict chromosome-specific regions
    sys.stderr.write('----region_dnacopy.R----\n')
    if not os.path.isfile(fnames['reg_tsv']):
        rd_command = [exec_path + '/region_dnacopy.R', fnames['pos_bed'], conf['sizes_file'], '20', '20']
        sys.stderr.write(' '.join(rd_command)+'\n')
        if not args.dry_run:
            p = subprocess.Popen(rd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
            (out, err) = p.communicate()
            if p.returncode != 0:
                raise Exception("region_dnacopy.R failed:\n\n%s" % (err))
    else:
        sys.stderr.write('%s chromosome region prediction file exists. OK!\n' % fnames['reg_tsv'])
    sys.stderr.write('----Complete!----\n')

    # Calculate sample stats if these do not exist. Less detailed than control_stats.
    sys.stderr.write('----sample_stats.py----\n')
    if not os.path.isfile(fnames['stat_file']):
        sys.stderr.write('Writing stats to %s\n' % fnames['stat_file'])
        ss_args = argparse.Namespace(sample=conf['sample'],genome=fnames['target_name'],
                            f_reads=conf['fastq_F_file'], r_reads=conf['fastq_R_file'], 
                            target_regions=None, output=None, dry_run = args.dry_run)
        sample_stats.main(ss_args)
    else:
        sys.stderr.write('%s statistics file exist. OK!\n' % fnames['stat_file'])
    sys.stderr.write('----Complete!----\n')
