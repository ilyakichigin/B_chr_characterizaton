#!/usr/bin/env python

import os
import sys
import argparse

from dopseq.tools import utils

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Remove primers and Illumina adapter from reads 
                    """
                    )
    parser.add_argument('sample', 
                        help="Sample name used as output prefix")
    parser.add_argument('fastq_F_file', 
                        help='Fastq file with forward reads')
    parser.add_argument('fastq_R_file', nargs='?', default=None, 
                        help='Optional: fastq file with reverse reads')
    parser.add_argument('-o', default=None, 
                        help='Output trimmed forward reads')
    parser.add_argument('-p', default=None, 
                        help='Output trimmed reverse reads')
    parser.add_argument('-l', default=None, 
                        help='Log file')
    parser.add_argument('-p', default=None, 
                        help='Output trimmed reverse reads')
    parser.add_argument('-d', '--dry_run', action='store_true', default=False, 
                        help='Print out commands and exit')
    parser.add_argument('--cutadapt_path', default='cutadapt', 
                        help='Path to cutadapt binary. Default: cutadapt')
    parser.add_argument('--trim_illumina', action='store_true', default=True, 
                        help='Trim illumina adapters. Default: True')
    parser.add_argument('--ampl', default='dop', 
                        help='Amplification protocol - used to remove specific'
                        ' primers. Valid values: dop, wga, none. '
                        'Default: dop.')
    parser.add_argument('--params', default='--trim-n --minimum-length 20', 
                        help='Additional parameters for cutadapt passed as '
                        'quoted string. '
                        'Default - trim terminal Ns and discard read pairs '
                        'with at least one read shorter than 20 (bowtie2 '
                        'seed length. For WGA it is sometimes useful to '
                        'increase error toleance, e.g. "-e 0.2"')

    return parser.parse_args()

def main(args):

    ca_reads = [args.o, args.p]
    if not args.o:
        ca_reads[0] = '%s.ca.R1.fastq.gz' % args.sample
    if not args.p:
        ca_reads[1] = '%s.ca.R1.fastq.gz' % args.sample 
    log_file = (args.l if args.l else args.sample + '.ca.log')

    # custom parameters
    cutadapt_opts = args.cutadapt_path + ' ' + args.params.strip("\'\"")
    # remove Illumna TruSeq adapters
    if args.trim_illumina:
        cutadapt_opts += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
        if args.fastq_R_file:
           cutadapt_opts += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    # remove primers from both ends of reads depending on the protocol
    if args.ampl == 'dop':
        cutadapt_opts += ' -n 3 -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG'
        if args.fastq_R_file:
            cutadapt_opts += ' -A CCACATNNNNNNCTCGAGTCGG -G CCGACTCGAGNNNNNNATGTGG'
    elif args.ampl == 'wga':
        cutadapt_opts += ' -n 3 -a CCAAACACACCCAACACAA -g TTGTGTTGGGTGTGTTTGG'
        if args.fastq_R_file:
            cutadapt_opts += ' -A CCAAACACACCCAACACAA -G TTGTGTTGGGTGTGTTTGG'
    elif args.ampl == 'none':
        pass #cutadapt_opts += ' '
    else:
         raise ValueError('Unknown amplification protocol %s. ' 
                          'Use dop, wga, or none instead.')

    # inputs and outputs. Reads left unpaired after trimming are discarded.
    if args.fastq_R_file:
        cutadapt_opts += ' -o %s -p %s %s %s' % (ca_reads[0], ca_reads[1], 
                                                 args.fastq_F_file, 
                                                 args.fastq_R_file)
    else:
        cutadapt_opts += ' -o %s %s' % (ca_reads[0], args.fastq_F_file)

    utils.run_command(cutadapt_opts, 
                      verbose=True, 
                      dry_run=args.dry_run,
                      outfile=log_file)

if __name__ == '__main__':
    main(parse_command_line_arguments())

