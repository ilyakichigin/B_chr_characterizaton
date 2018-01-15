#!/usr/bin/env python

import os
import sys
import argparse

from dopseq.tools import utils

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Align fastq files to reference genome with bowtie2 or bwa
                    Required program: bowtie2 (tested on v.2.1.0, 2.2.4). 
                    """
                    )

    parser.add_argument('fastq_F_file', help='Fastq file with forward reads')
    parser.add_argument('fastq_R_file', nargs='?', default=None,
                        help='Optional: fastq file with reverse reads')
    parser.add_argument('-r', '--reference_genome', 
                        help='Path to reference genome unpacked fasta')
    parser.add_argument('-o', '--out_bam', default='test.bam',
                        help='Output BAM file name')
    parser.add_argument('-a', '--aligner', default='bwa',
                        help='Aligner to use. Possible values: '
                             'bt2 (bowtie2), bwa.')
    parser.add_argument('-p', '--aligner_path', default=None, 
                        help='Path to aligner program')
    parser.add_argument('-b', '--aligner_args', default='mem -t 1', 
                        help='Additional parameters for aligner specified as '
                        'quoted string. For reference, see alignmer manuals.')
    parser.add_argument('-d', '--dry_run', action='store_true', default=False, 
                        help='Print out commands and exit')

    return parser.parse_args()


def main(args):    

    # outputs
    assert args.out_bam.endswith('.bam')
    out_sam = '%s.sam' % (args.out_bam[:-4]) 
    log_file = '%s.%s.log' % (args.out_bam[:-4], args.aligner) 
    
    if args.aligner == 'bt2':
        if not args.aligner_path:
            args.aligner_path = 'bowtie2'
        # index
        if not os.path.isfile(args.reference_genome + '.1.bt2'):
            bt2_index_command = '%s-build %s %s' % (args.aligner_path, 
                                                    args.reference_genome, 
                                                    args.reference_genome)
            utils.run_command(bt2_index_command, dry_run=args.dry_run, 
                              outfile=None, errfile=None)
        # align
        if not args.fastq_R_file:
            bt2_align_command = '%s %s -U %s -x %s -S %s' % (
                args.aligner_path, args.aligner_args, 
                args.fastq_F_file, args.reference_genome, out_sam)
        else:
            bt2_align_command = '%s %s -1 %s -2 %s -x %s -S %s' % (
                args.aligner_path, args.aligner_args, 
                args.fastq_F_file, args.fastq_R_file, 
                args.reference_genome, out_sam)
        utils.run_command(bt2_align_command, dry_run=args.dry_run, 
                              outfile=None, errfile=log_file)

    elif args.aligner == 'bwa':
        if not args.aligner_path:
            args.aligner_path = 'bwa'
        if not os.path.isfile(args.reference_genome + '.bwt'):
            bwa_index_command = '%s index %s' % (args.aligner_path, 
                                                 args.reference_genome)
            utils.run_command(bwa_index_command, dry_run=args.dry_run, 
                              outfile=None, errfile=None)
        if not args.fastq_R_file:
            bwa_align_command = '%s %s %s %s' % (args.aligner_path, 
                args.aligner_args, args.reference_genome, args.fastq_F_file)
        else:
            bwa_align_command = '%s %s %s %s %s' % ( 
                args.aligner_path, args.aligner_args, args.reference_genome, 
                args.fastq_F_file, args.fastq_R_file)
        utils.run_command(bwa_align_command, dry_run=args.dry_run, 
                              outfile=out_sam, errfile=log_file)

    elif args.aligner == 'bbm':
        raise ValueError('bbmap alignment support discontinued.'
                         'Please use bowtie2 or bwa instead.')

    else:
        raise ValueError('Invalid aligner!')

    utils.bam_sort(out_sam, args.out_bam, dry_run=args.dry_run)
    if not args.dry_run:
        os.unlink(out_sam)

if __name__ == '__main__':
    main(parse_command_line_arguments())
