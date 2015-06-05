#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Pipeline for processing of sequencing data of DOP-PCR libraries from isolated chromosomes.
                    See config for process description, inputs and outputs. 
                    """
                    )
    parser.add_argument("config_file", help="input configuration file")
    
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

    # Does not stop on errors from within
    sys.stderr.write(' '.join(command)+'\n') 
    if run:
        process = subprocess.Popen(command) 
        process.wait()

def run_script_to_file(command, outfile, run=False):
    
    sys.stderr.write(' '.join(command)+'\n') 
    if run:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out,err) = process.communicate()
        sys.stderr.write(err)
        with open(outfile, 'w') as of:
            of.write(out) 
        process.wait()
  
    
if __name__ == '__main__':
    args = parse_command_line_arguments()
    conf = parse_config(args.config_file)
    sample = conf['sample']
    assert len(sample) > 0 # sample name not empty
    path_to_exec = conf['path_to_exec']
    target_name = conf["target_genome"].split('/')[-1]
    contam_name = conf["contam_genome"].split('/')[-1]

    if not path_to_exec.endswith('/'):
        path_to_exec += '/'
       
    if conf['do_fastq_to_bam'] == 'True':
        command = [path_to_exec + "fastq_to_bam.py","-t",conf["target_genome"],"-c",conf["contam_genome"],
                    "--path_to_cutadapt",conf["path_to_cutadapt"],"--path_to_bowtie2",conf["path_to_bowtie2"],
                    "-p",conf["proc_bowtie2"],"--wga",conf["wga"],conf['fastq_F_file'],conf['fastq_R_file'],sample]
        run_script(command, run=True)
    if conf['do_contam_filter'] == 'True':
        target_sam = '.'.join([sample,target_name,'sam'])
        contam_sam = '.'.join([sample,contam_name,'sam'])
        command = [path_to_exec + "contam_filter.py","-m",conf["min_quality"],"-a",
                    target_sam,contam_sam]
        run_script(command, run=True)
    if conf['do_bam_to_beds'] == 'True':
        in_bam = '.'.join([sample,target_name,'filter','bam'])
        command = [path_to_exec + "bam_to_beds.py","--path_to_bedtools",conf["path_to_bedtools"],in_bam]
        run_script(command, run=True)
    if conf['do_region_dnacopy'] == 'True':
        pos_bed = '.'.join([sample,target_name,'filter','pos','bed'])
        command = [path_to_exec + "region_dnacopy.R",pos_bed,conf["sizes_file"]]
        run_script(command, run=True)
