#!/usr/bin/env python

import subprocess
import argparse
import os.path
import pysam

import sys
exec_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"..","exec"))
sys.path.append(exec_path)

import fastq_clean
import fastq_clean_se
import fastq_to_bam
import contam_filter
import bam_to_beds
import control_stats
import sample_stats

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Pipeline for processing of sequencing data of DOP-PCR libraries from isolated chromosomes.
                    See config for process description, inputs and outputs. 
                    """
                    )
    parser.add_argument("config_file", help="input configuration file")
    #parser.add_argument("-d", "--dry-run", action="store_true",
    #                    help="Check all dependencies and print out all commands")
    
    return parser.parse_args()

def parse_config(config_file):
    
    conf = dict()
    
    with open(config_file) as infile:
        for line in infile:
            line = line[:line.find('#')] # removes comments
            line_list = line.split('=')
            if len(line_list) == 2:            
                conf[line_list[0]]=line_list[1].strip(' ')
            #elif len(line_list) == 1:
            #    raise Exception(line+"\n This line in config has no parameter setting and is not comment")
            elif len(line_list) > 2:
                raise Exception(line+"\n This line in config has too many '=' symbols")

    return conf

def run_script(command, run=False):

    # Err: Does not stop on errors from within
    sys.stderr.write(' '.join(command)+'\n') 
    process = subprocess.Popen(command) 
    process.wait()
    if process.returncode != 0: # error raised
        sys.exit(1)

if __name__ == '__main__':
    args = parse_command_line_arguments()
    conf = parse_config(args.config_file)
    #dry_run = args.dry_run
    # !need to add executables check!
    f_trim_fq = conf['sample']+'.ca.R1.fastq'
    if "fastq_R_file" in conf.keys():
        r_trim_fq = conf['sample']+'.ca.R2.fastq'
    target_name = conf["target_genome"].split('/')[-1]
    target_name = target_name.split('.')[0] + '_' + conf["aln"]
    base_name = '.'.join([conf['sample'],target_name,'filter'])
    filtered_bam_file = base_name + '.bam'
    if not os.path.isfile(filtered_bam_file): 

        # Step 1. fastq_clean if trimmed read fastq do not exists
        # paired-end
        if not os.path.isfile(f_trim_fq) and 'fastq_R_file' in conf.keys():
            fc_args = argparse.Namespace(fastq_F_file=conf['fastq_F_file'],fastq_R_file=conf['fastq_R_file'],
                                         sample_name=conf['sample'],path_to_cutadapt='cutadapt',
                                         ampl=conf["ampl"],params=conf["cutadapt_args"],
                                         delimiter=' ', rename_only=False)
            assert os.path.isfile(conf['fastq_F_file'])
            assert os.path.isfile(conf['fastq_R_file'])
            sys.stderr.write('----fastq_clean.py----\n')
            fastq_clean.main(fc_args)
            sys.stderr.write('----Complete!----\n')
        # single-end
        if not os.path.isfile(f_trim_fq) and 'fastq_R_file' not in conf.keys():
            fcse_args = argparse.Namespace(fastq_file=conf['fastq_F_file'],
                                         sample_name=conf['sample'],path_to_cutadapt='cutadapt',
                                         ampl=conf["ampl"],params=conf["cutadapt_args"])
            assert os.path.isfile(conf['fastq_F_file'])
            sys.stderr.write('----fastq_clean_se.py----\n')
            fastq_clean_se.main(fcse_args)
            sys.stderr.write('----Complete!----\n')

        # Step 2a. Align to target genome
        target_sam_file = '.'.join([conf['sample'],target_name,'sam'])
        if conf['aln'] == 'bt2':
            path_to_aligner = conf['path_to_bowtie2']
            aligner_args = conf['bowtie2_args']
        elif conf['aln'] == 'bbm':
            path_to_aligner = conf['path_to_bbmap']
            aligner_args = conf['bbmap_args']
        else:
            raise Exception('Invalid aligner!')
        fb2t_args = argparse.Namespace(sample=conf['sample'], aligner=conf["aln"], reference_genome=conf["target_genome"],
                                 path_to_aligner=path_to_aligner[1:-1], aligner_args=aligner_args[1:-1])
        sys.stderr.write('----fastq_to_bam.py----\n')
        fastq_to_bam.main(fb2t_args)

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
            cf_args = argparse.Namespace(target_file=target_sam_file,contam_file=contam_sam_file,
                                         min_quality=20)
            sys.stderr.write('----contam_filter.py----\n')
            contam_filter.main(cf_args)
            sys.stderr.write('----Complete!----\n')
        else: # no contam filter
            sys.stderr.write('Sorting and indexing %s, resulting file %s\n' % (target_sam_file, filtered_bam_file))
            sys.stderr.write('samtools view -bSq 20 %s | samtools sort -n - > %s\n'%(target_sam_file, filtered_bam_file))
            p1 = subprocess.Popen(('samtools', 'view', '-bSq', '20', target_sam_file), stdout=subprocess.PIPE)
            #with p1.stdout, open(filtered_bam_file, 'w') as outfile:
            p2 = subprocess.Popen(('samtools', 'sort', '-', filtered_bam_file[:-4]), stdin=p1.stdout, stdout=subprocess.PIPE)
            status=[p1.wait(),p2.wait()]
            if p1.returncode != 0 or p1.returncode != 0:
                sys.exit(1)
            pysam.index(filtered_bam_file)
            #os.remove(target_sam_file)
            sys.stderr.write('----Complete!----\n')

    # Step 4. Convert bam_to_beds with reads and positions if these files do not exist.
    reads_bed_file = base_name+'.reads.bed'
    pos_bed_file = base_name+'.pos.bed'
    if (not os.path.isfile(reads_bed_file) and not os.path.isfile(pos_bed_file)):
        btb_args = argparse.Namespace(bam_file=filtered_bam_file,path_to_bedtools='bedtools')
        sys.stderr.write('----bam_to_beds.py----\n')
        #try:
        bam_to_beds.main(btb_args)
        #except:
        #    sys.exit(1)
        sys.stderr.write('----Complete!----\n')
    # Step 4a. Calculate control_stats if these do not exist.
    stat_file = base_name+'.chrom.tsv'
    if not os.path.isfile(stat_file):
        cs_args = argparse.Namespace(bed_basename='.'.join([conf['sample'],target_name,'filter']),
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
    control_plot_file = base_name+'.chrom.pdf'
    if not os.path.isfile(control_plot_file):
        cp_command = [exec_path+'/control_plots.R',pos_bed_file,conf['sizes_file']]
        sys.stderr.write('----control_plots.R----\n')
        run_script(cp_command)
        sys.stderr.write('----Complete!----\n')
    # Step 5. Perform region_dnacopy.R if regions do not exist.
    regions_table = base_name+'.reg.tsv'
    regions_plot = base_name+'.reg.pdf'
    if not os.path.isfile(regions_table) or not os.path.isfile(regions_plot):
        # pdf height and width increased to get readable plot for all chromosomes 
        rd_command = [exec_path+'/region_dnacopy.R',pos_bed_file,conf['sizes_file'],'20','20']
        sys.stderr.write('----region_dnacopy.R----\n')
        run_script(rd_command)
        sys.stderr.write('----Complete!----\n')
    '''
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