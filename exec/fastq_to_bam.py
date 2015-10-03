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

    parser.add_argument("sample", help="sample name - prefix of input files processed with cutadapt")

    parser.add_argument("-t", "--target_genome", help="base name of reference genome bowtie2 index (e.g. canFam3, bosTau7)")

    parser.add_argument("-c", "--contam_genome", default="hg19", help="base name for contamination genome bowtie2 index (default hg19)")

    parser.add_argument("--path_to_bowtie2", default="bowtie2", help="path to bowtie2 binary")

    parser.add_argument("-b", "--bowtie2_args", default="-p 1", help="Additional parameters for bowtie2 specified as quoted string. For reference, see bowtie2 manual.")

    return parser.parse_args()

def run_bt2(command, log_file):
    
    # run command logging to a file
    sys.stderr.write(command+'\n') 
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    (out, err) = process.communicate()
    with open(log_file, 'w') as log:
        log.write(out)
    # ignore bowtie2 warnings for too short reads. All sterr is stored in memory!         
    #for line in err.splitlines(True):
    #    if ('Warning: skipping mate' not in line) and ('Warning: minimum score function' not in line):
    #        
    sys.stderr.write(err)
    if process.returncode != 0:
        sys.exit()

def main(args):    

    # assign input filenames
    f_ca_fq = args.sample + '.ca.R1.fastq'
    r_ca_fq = args.sample + '.ca.R2.fastq'
    # assign output filenames
    target = args.target_genome.split('/')[-1]
    contam = args.contam_genome.split('/')[-1]
    target_sam = args.sample + '.' + target + '.sam' # simplified name: sample.genome.sam. 'ca' and 'pe not included
    contam_sam = args.sample + '.' + contam + '.sam'
    
    # alignment to genomes - if not already done
    command_base = args.path_to_bowtie2 + ' ' + args.bowtie2_args + ' -1 ' + f_ca_fq + ' -2 ' + r_ca_fq 
    if (not os.path.isfile(target_sam)) or (os.path.getsize(target_sam) == 0):
        command = command_base + ' -x %s -S %s' % (args.target_genome, target_sam)
        log_file = target_sam[:-3]+'bt2.log'
        run_bt2(command, log_file)
    else:
        print 'Alignment to target genome exists. OK!'
    if (not os.path.isfile(contam_sam)) or (os.path.getsize(contam_sam) == 0):
        command = command_base + ' -x %s -S %s' % (args.contam_genome, contam_sam)
        log_file = contam_sam[:-3]+'bt2.log'
        run_bt2(command, log_file)
    else:
        print 'Alignment to contamination genome exists. OK!'

if __name__ == '__main__':
    main(parse_command_line_arguments())
