#!/usr/bin/env python

import subprocess
import sys
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Aligns pair of fastq read files to target and contamination genomes with preliminary DOP primer and Illumina adapter removal, prints out programs used to stdout. 
                    Required programs: cutadapt (tested on v.1.6), bowtie2 (tested on v.2.1.0, 2.2.4). 
                    1) Change read names to include '1' for forward and '2' for reverse read (Illumina-specific).
                    2) Cut Illumina adapters and DOP primers.
                    3) Perform paired end mapping to reference and contamination genomes.
                    """
                    )
    parser.add_argument("fastq_F_file", help="fastq file with forward reads (.fastq)")

    parser.add_argument("fastq_R_file", help="fastq file with reverse reads (.fastq)")

    parser.add_argument("sample_name", help="output name of sample (can be anything, in our case CFA12, VVUB and so on)")

    parser.add_argument("-t", "--target_genome", help="base name of reference genome (e.g. canFam3, bosTau7)")

    parser.add_argument("-c", "--contam_genome", default="hg19", help="base name for contamination genome (default hg19)")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--path_to_bowtie2", default="bowtie2", help="path to bowtie2 binary")

    parser.add_argument("-p", "--proc_bowtie2", default="1", help="number of processors allocated for bowtie2. Default - 1.")

    return parser.parse_args()
   
def rename_reads(in_file_name):

    # Rename reads in fastq files - remove space between cluster name and number: 1 for F, 2 for R. 

    with open(in_file_name, 'rU') as in_file:
        in_file_base = in_file_name.split('/')[-1] # output to current folder
        out_file_name = '.'.join(in_file_base.split('.')[:-1]+['rn','fastq'])
        sys.stdout.write("Renaming reads in %s. Writing to %s.\n" % (in_file_name,out_file_name))
        with open(out_file_name, 'w') as out_file:
            for line in in_file:
                line_list = line.split(' ')
                if line.startswith('@') and len(line_list) > 1:
                    out_line = line_list[0] + line_list[1] # remove separating space
                else:
                    out_line = line
                out_file.write(out_line)

    return out_file_name    

if __name__ == '__main__':
    args = parse_command_line_arguments()

    assert args.fastq_F_file.endswith('.fastq') and args.fastq_R_file.endswith('.fastq') # does not accept gzipped and improperly named files

    # rename reads
    forward_rn_fq = rename_reads(args.fastq_F_file) # Filename returned. rn = renamed
    reverse_rn_fq = rename_reads(args.fastq_R_file)
    
    # assign output filenames
    f_ca_fq_name = args.sample_name + '.F.ca.fastq'
    r_ca_fq_name = args.sample_name + '.R.ca.fastq'
    target_sam_name = args.sample_name + '.' + args.target_genome.split('/')[-1] + '.sam' # simplified name: sample.genome.sam. 'ca' and 'pe not included (compared to previous versions)
    contam_sam_name = args.sample_name + '.' + args.contam_genome.split('/')[-1] + '.sam'
    
    # generate commands
    # * supress bowtie2 warnings - too many for short reads
    command_list = [
            (args.path_to_cutadapt + ' -a AGATCGGAAGAGC -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3 -o ' + f_ca_fq_name + ' ' + forward_rn_fq),
            (args.path_to_cutadapt + ' -a AGATCGGAAGAGC -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3 -o ' + r_ca_fq_name + ' ' + reverse_rn_fq),
            (args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x ' + args.target_genome + ' -1 ' + f_ca_fq_name + ' -2 ' + r_ca_fq_name + ' -S ' + target_sam_name),
            (args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x ' + args.contam_genome + ' -1 ' + f_ca_fq_name + ' -2 ' + r_ca_fq_name + ' -S ' + contam_sam_name)
            ]
    
    # run
    for command in command_list:
        sys.stdout.write(command+'\n') 
        process = subprocess.Popen(command.split()) 
        '''        
        # ignore bowtie2 warnings of too short reads. All sterr is stored in memory!         
        err = process.communicate()[1]
        for line in err:
            if ('Warning: skipping mate' not in line) and ('Warning: minimum score function' not in line):
                sys.stderr.write(line)
        '''        
        process.wait()
    
    sys.stdout.write("Complete!\n")
