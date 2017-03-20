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
#import control_stats
#import sample_stats

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
    names['f_trim_fq'] = sample + '.ca.R1.fastq'
    names['r_trim_fq'] = sample + '.ca.R2.fastq'
    names['target_path'] = target_path
    names['contam_path'] = contam_path
    names['target_name'] = target_path.split('/')[-1].split('.')[0]
    names['target_index_prefix'] = '/'.join(target_path.split('/')[:-1]) + names['target_name']
    names['contam_name'] = contam_path.split('/')[-1].split('.')[0]
    names['contam_index_prefix'] = '/'.join(target_path.split('/')[:-1]) + names['contam_name']
    names['target_sam'] = '.'.join([sample,names['target_name'],'sam'])
    names['contam_sam'] = '.'.join([sample,names['contam_name'],'sam'])
    names['filtered_bam'] = '.'.join([sample,names['target_name'],'filter','bam'])
    names['pos_bed'] = '.'.join([sample,names['target_name'],'pos','bed'])
    names['reg_tsv'] = '.'.join([sample,names['target_name'],'reg','tsv'])
    names['reg_pdf'] = '.'.join([sample,names['target_name'],'reg','pdf'])

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
        elif conf['aln'] == 'bbm':
            aligner = conf['bbmap_path']
            aligner_args = conf['bbmap_args']
        else:
            raise Exception('Invalid aligner!')
        for g in ('target', 'contam'):
            fa_args = argparse.Namespace(fastq_F_file=fnames['f_trim_fq'], 
                fastq_R_file=fnames['r_trim_fq'], aligner=conf['aln'], 
                reference_genome=fnames[g + '_path'], path_to_aligner=aligner_path, 
                aligner_args=aligner_args, dry_run = args.dry_run)
            fastq_aln.main(fa_args)
        sys.stderr.write('----Complete!----\n')

        sys.stderr.write('----contam_filter.py----\n') # filter contamination
        cf_args = argparse.Namespace(target_sam=fnames['target_sam'], contam_sam=fnames['contam_sam'],
                                     min_quality=conf['min_mapq'], dry_run = args.dry_run)
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
'''
    base_name = '.'.join([conf['sample'], target_name,'filter'])
    filtered_bam_file = base_name + '.bam'
    if not os.path.isfile(filtered_bam_file): 

        # Step 1. fastq_clean if trimmed read fastq do not exists
        # paired-end


        

        if "contam_genome" in conf.keys():  # contam filter
            contam_name = conf["contam_genome"].split('/')[-1]
            contam_name = contam_name.split('.')[0] + '_' + conf["aln"]
            contam_sam_file = '.'.join([conf['sample'],contam_name,'sam'])
            # Step 2b. Align to contamination genome
            fb2c_args = argparse.Namespace(sample=conf['sample'],  aligner=conf["aln"], reference_genome=conf["contam_genome"], 
                                         path_to_aligner=path_to_aligner[1:-1], aligner_args=aligner_args[1:-1])
            fastq_to_bam.main(fb2c_args)
            sys.stderr.write('----Complete!----\n')

            # Step 3. contam_filter - remove contamination from the specified genome
            cf_args = argparse.Namespace(target_file=target_sam_file, contam_file=contam_sam_file,
                                         min_quality=conf['min_mapq'])
            sys.stderr.write('----contam_filter.py----\n')
            contam_filter.main(cf_args)
            sys.stderr.write('----Complete!----\n')
        else: # no contam filter
            sys.stderr.write('Sorting and indexing %s, resulting file %s\n' % (target_sam_file, filtered_bam_file))
            view_command = 'samtools view -bSq %s %s' % (conf['min_mapq'], target_sam_file)
            sort_command = 'samtools sort - %s' % (filtered_bam_file[:-4])
            sys.stderr.write('%s | %s\n' % (view_command, sort_command))
            p1 = subprocess.Popen(view_command.split(), stdout=subprocess.PIPE)
            p2 = subprocess.Popen(sort_command.split(), stdin=p1.stdout, stdout=subprocess.PIPE)
            status = [p1.wait(),p2.wait()]
            if p1.returncode != 0 or p1.returncode != 0:
                sys.exit(1)
            pysam.index(filtered_bam_file)
            #os.remove(target_sam_file)
            sys.stderr.write('----Complete!----\n')

    # Step 4. Convert bam_to_beds with reads and positions if these files do not exist.
    reads_bed_file = base_name + '.reads.bed'
    pos_bed_file = base_name + '.pos.bed'
    if (not os.path.isfile(reads_bed_file) and not os.path.isfile(pos_bed_file)):
        btb_args = argparse.Namespace(bam_file=filtered_bam_file, path_to_bedtools='bedtools')
        sys.stderr.write('----bam_to_beds.py----\n')
        #try:
        bam_to_beds.main(btb_args)
        #except:
        #    sys.exit(1)
        sys.stderr.write('----Complete!----\n')
    # Step 4a. Calculate control_stats if these do not exist.
    stat_file = base_name + '.chrom.tsv'
    if not os.path.isfile(stat_file):
        cs_args = argparse.Namespace(bed_basename='.'.join([conf['sample'], target_name,'filter']),
                                     sizes_file=conf['sizes_file'])
        # redirect stdout to file and then return it back
        sys.stderr.write('----control_stats.py----\n')        
        saveout = sys.stdout        
        with open(stat_file,'w') as f:
            sys.stdout = f # replace stdout with file
            #try:
            control_stats.main(cs_args)
            #except:
            #    sys.exit(1)
            sys.stdout = saveout # return std
        sys.stderr.write('----Complete!----\n')
    # Step 4b. Draw control_plots.R if these do not exist.
    control_plot_file = base_name + '.chrom.pdf'
    if not os.path.isfile(control_plot_file):
        cp_command = [exec_path + '/control_plots.R', pos_bed_file,conf['sizes_file']]
        sys.stderr.write('----control_plots.R----\n')
        run_script(cp_command)
        sys.stderr.write('----Complete!----\n')
    # Step 5. Perform region_dnacopy.R if regions do not exist.
    regions_table = base_name + '.reg.tsv'
    regions_plot = base_name + '.reg.pdf'
    if not os.path.isfile(regions_table) or not os.path.isfile(regions_plot):
        # pdf height and width increased to get readable plot for all chromosomes 
        rd_command = [exec_path + '/region_dnacopy.R', pos_bed_file,conf['sizes_file'], '20', '20']
        sys.stderr.write('----region_dnacopy.R----\n')
        run_script(rd_command)
        sys.stderr.write('----Complete!----\n')
    
    # Step 6. Calculate sample stats if these do not exist. Less detailed than control_stats.
    stat_file = conf['sample'] + '.stats.txt'
    if not os.path.isfile(stat_file):
        
        if "contam_genome" in conf.keys():
            contam_name = conf["contam_genome"].split('/')[-1]
            ss_args = argparse.Namespace(sample=conf['sample'],
                            F_reads=conf['fastq_F_file'], R_reads=conf['fastq_R_file'],
                            t_genome=target_name,c_genome=contam_name)
        else:
            ss_args = argparse.Namespace(sample=conf['sample'],
                            F_reads=conf['fastq_F_file'], R_reads=conf['fastq_R_file'],
                            t_genome=target_name,c_genome='')
        sys.stderr.write('----sample_stats.py----\n')
        sample_stats.main(ss_args)
        #try:
            
        #except:
        #    sys.exit(1)
        sys.stderr.write('----Complete!----\n')
    '''