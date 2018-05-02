#!/usr/bin/env python

import sys
from gzip import GzipFile
import argparse

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Convert paired reads from Illumina 1.8+ fastq processed 
                    with cutadapt to interleaved fasta suitable as input for 
                    RepeatExplorer with 'All sequence reads are paired' option.
                    """
                    )
    parser.add_argument('f_fastq', # args.cutadapt_prefix?
                        help='file with forward reads')
    parser.add_argument('r_fastq', 
                        help='file with forward reads')
    parser.add_argument('-o', '--out_fasta', default=None,
                        help='output fasta filename')
    parser.add_argument('-r', '--rename', action='store_true',
                        help='rename reads to consecutive numbers')
    parser.add_argument('-d', '--dry_run', action='store_true',
                        help='print output file name and exit')
    return parser.parse_args()
   
def main(args):
    """From two fastq files generate interleaved fasta file with added '/1' 
    for forward and '/2' for reverse read
    """
    
    if not args.out_fasta:
        args.out_fasta = args.cutadapt_prefix + '.ca.re.fasta' # args.cutadapt_prefix uninitialized! 
    rename = args.rename
    
    # minimum read length: word length in megablast used for cluster annotation
    min_read_length = 18 
    i = 0 # read number
    k = 0 # line number
    f_name = '' # forward read name
    r_name = '' # reverse read name
    sys.stderr.write('Output file: %s\n' % args.out_fasta)
    if args.dry_run:
        return 0

    with GzipFile(args.f_fastq, 'r') as f_file, \
         GzipFile(args.r_fastq, 'r') as r_file, \
         open(args.out_fasta, 'w') as out:
        for f_line in f_file:
            r_line = r_file.next()
            # read name
            if f_line.startswith('@') and (k % 4) == 0: # get read name
                i += 1
                if rename:
                    f_name = '>%s_%d/1\n' % (args.cutadapt_prefix, i)
                    r_name = '>%s_%d/2\n' % (args.cutadapt_prefix, i)
                else:
                    f_name = '>%s/1\n' % f_line[1:-1]
                    r_name = '>%s/2\n' % r_line[1:-1]
            # sequence 
            elif f_name != '' and r_name != '':
                if (len(f_line) > min_read_length + 1 
                        and len(r_line) > min_read_length + 1): 
                    out.write(f_name + f_line + r_name + r_line)
                f_name = ''
                r_name = ''
            
            k += 1

if __name__ == '__main__':
    main(parse_command_line_arguments())

    
