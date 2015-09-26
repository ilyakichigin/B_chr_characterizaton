#!/usr/bin/env python

import subprocess
import sys
import os
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes primers and Illumina adapter, prints out programs used to stdout. 
                    Required program: cutadapt (tested on v.1.6)
                    1) Change read names to include '1' for forward and '2' for reverse read (Illumina-specific).
                    2) Cut Illumina adapters and DOP or WGA primers.
                    3) Perform paired end mapping to reference and contamination genomes.
                    """
                    )
    parser.add_argument("fastq_F_file", help="fastq file with forward reads (.fastq)")

    parser.add_argument("fastq_R_file", help="fastq file with reverse reads (.fastq)")

    parser.add_argument("sample_name", help="sample name - used as output prefix")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--ampl", default="dop", help="Amplification protocol - used to remove specific primers. Possible values: dop, wga, none")

    return parser.parse_args()
   
def rename_reads(in_file_name):

    # Rename reads in fastq files - remove space between cluster name and number: 1 for F, 2 for R. 

    with open(in_file_name, 'rU') as in_file:
        in_file_base = in_file_name.split('/')[-1] # output to current folder
        out_file_name = '.'.join(in_file_base.split('.')[:-1]+['rn','fastq'])
        sys.stderr.write("Renaming reads in %s. Writing to %s.\n" % (in_file_name,out_file_name))
        with open(out_file_name, 'w') as out_file:
            for line in in_file:
                line_list = line.split(' ')
                if line.startswith('@') and len(line_list) > 1:
                    out_line = line_list[0] + line_list[1] # remove separating space
                else:
                    out_line = line
                out_file.write(out_line)

    return out_file_name

def main(args):
    
    assert args.fastq_F_file.endswith('.fastq') and args.fastq_R_file.endswith('.fastq') # do not accept gzipped and improperly named files
    # assign output filenames
    f_ca_fq_name = args.sample_name + '.F.ca.fastq'
    r_ca_fq_name = args.sample_name + '.R.ca.fastq'

    # rename and cutadapt if not already done
    if (not os.path.isfile(f_ca_fq_name) and not os.path.isfile(r_ca_fq_name)):    
        # rename reads
        forward_rn_fq = rename_reads(args.fastq_F_file) # Filename returned. rn = renamed
        reverse_rn_fq = rename_reads(args.fastq_R_file)
        # remove Illumina unversal adapter and primers depending on the protocol
        if args.ampl == 'dop':
            cutadapt_opts = ' -a AGATCGGAAGAGC -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3'
        elif args.ampl == 'wga':
            cutadapt_opts = ' -a AGATCGGAAGAGC -g TTGTGTTGGGTGTGTTTGG -a CCAAACACACCCAACACAA -n 3'
        elif args.ampl == 'none':
            cutadapt_opts = ' -a AGATCGGAAGAGC'
        else:
             raise Exception('Unknown amplification protocol. Known ones - dop, wga, none')
        command_list = [
                (args.path_to_cutadapt + cutadapt_opts + ' -o ' + f_ca_fq_name + ' ' + forward_rn_fq),
                (args.path_to_cutadapt + cutadapt_opts + ' -o ' + r_ca_fq_name + ' ' + reverse_rn_fq),
                ('rm ' + forward_rn_fq + ' ' + reverse_rn_fq)]
        for command in command_list:
            sys.stderr.write(command+'\n') 
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
            (out, err) = process.communicate()
            sys.stdout.write(out)
    else:
        command_list = []
        print 'Files with removed adapters exist. OK!'


if __name__ == '__main__':
    main(parse_command_line_arguments())

