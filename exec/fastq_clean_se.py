#!/usr/bin/env python

import subprocess
import sys
import os
import argparse
from gzip import GzipFile


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes primers and Illumina adapter in SE mode, prints out programs used to stderr. 
                    Required program: cutadapt (tested on v.1.8)
                    Output files: sample.ca.fastq, sample.ca.log. 
                    """
                    )
    parser.add_argument("fastq_file", help="fastq file with single-end reads (.fastq or .fastq.gz)")

    parser.add_argument("--sample_name", help="sample name - used as output prefix")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--ampl", default="dop", help="Amplification protocol - used to remove specific primers. Possible values: dop, wga, none")

    parser.add_argument("--params", default="--trim-n --minimum-length 20", help="Additional parameters for cutadapt. \
        Default - trim terminal Ns and discard read pairs with at least one read shorter than 20 (bowtie2 seed length). \
        For WGA it is sometimes useful to increase error toleance with -e")

    return parser.parse_args()

def main(args):

    ca_fq = args.sample_name + '.ca.fastq'
    log_file = args.sample_name + '.ca.log'
    if not os.path.isfile(ca_fq):
        
        # common options: remove Illumna TruSeq adapters, apply additional parameters.
        cutadapt_opts = args.path_to_cutadapt + ' ' + args.params.strip("\'\"") + ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
                                                
        # remove primers from both ends of reads depending on the protocol
        if args.ampl == 'dop':
            cutadapt_opts += ' -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3'
        elif args.ampl == 'wga':
            cutadapt_opts += ' -a CCAAACACACCCAACACAA -g TTGTGTTGGGTGTGTTTGG -n 3'
        elif args.ampl == 'none':
            cutadapt_opts += ' '
        else:
             raise Exception('Unknown amplification protocol. Known ones - dop, wga, none')
        # inputs and outputs.
        cutadapt_opts += ' -o %s %s' % (ca_fq, args.fastq_file)
        with open(log_file, 'w') as log:
            sys.stderr.write(cutadapt_opts+'\n') 
            process = subprocess.Popen(cutadapt_opts.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
            (out, err) = process.communicate()
            log.write(out)
            sys.stderr.write(err)
            if process.returncode != 0:
                sys.exit()
    else:
        print 'Files with removed adapters exist. OK!'


if __name__ == '__main__':
    main(parse_command_line_arguments())

