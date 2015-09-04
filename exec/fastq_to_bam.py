#!/usr/bin/env python

import subprocess
import sys
import os
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Aligns paired fastq files (preliminary processed with fastq_clean.py) to target and contamination genomes, prints out programs used to stdout. 
                    Required program: bowtie2 (tested on v.2.1.0, 2.2.4). 
                    """
                    )

    parser.add_argument("sample_name", help="sample name - prefix of input files processed with cutadapt")

    parser.add_argument("-t", "--target_genome", help="base name of reference genome bowtie2 index (e.g. canFam3, bosTau7)")

    parser.add_argument("-c", "--contam_genome", default="hg19", help="base name for contamination genome bowtie2 index (default hg19)")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--path_to_bowtie2", default="bowtie2", help="path to bowtie2 binary")

    parser.add_argument("-p", "--proc_bowtie2", default="1", help="number of processors allocated for bowtie2. Default - 1.")

    return parser.parse_args()

def main(args):    

    # assign input filenames
    f_ca_fq_name = args.sample_name + '.F.ca.fastq'
    r_ca_fq_name = args.sample_name + '.R.ca.fastq'
    # assign output filenames
    target_name = args.target_genome.split('/')[-1]
    contam_name = args.contam_genome.split('/')[-1]
    target_sam_name = args.sample_name + '.' + target_name + '.sam' # simplified name: sample.genome.sam. 'ca' and 'pe not included
    contam_sam_name = args.sample_name + '.' + contam_name + '.sam'
    
    # alignment to genomes - if not already done
    command_list = []
    if (not os.path.isfile(target_sam_name)) or (os.path.getsize(target_sam_name) == 0):
        command_list.append(args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x ' + args.target_genome + ' -1 ' + f_ca_fq_name + ' -2 ' + r_ca_fq_name + ' -S ' + target_sam_name)
    else:
        print 'Alignment to target genome exists. OK!'
    if (not os.path.isfile(contam_sam_name)) or (os.path.getsize(contam_sam_name) == 0):
        command_list.append(args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x ' + args.contam_genome + ' -1 ' + f_ca_fq_name + ' -2 ' + r_ca_fq_name + ' -S ' + contam_sam_name)
    else:
        print 'Alignment to contamination genome exists. OK!'
    # run
    for command in command_list:
        sys.stderr.write(command+'\n') 
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        (out, err) = process.communicate()
        sys.stdout.write(out)
        # ignore bowtie2 warnings of too short reads. All sterr is stored in memory!         
        for line in err.splitlines(True):
            if ('Warning: skipping mate' not in line) and ('Warning: minimum score function' not in line):
                sys.stderr.write(line)

if __name__ == '__main__':
    main(parse_command_line_arguments())
