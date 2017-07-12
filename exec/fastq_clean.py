#!/usr/bin/env python

import subprocess
import sys
import os
import argparse
from gzip import GzipFile


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes primers and Illumina adapter in PE mode, prints out programs used to stderr. 
                    Required program: cutadapt (tested on v.1.8)
                    1) Rename reads to include '/1' for forward and '/2' for reverse read (for cutadapt).
                    2) Trim Illumina TruSeq adapters, as well as DOP or WGA primers. Write log.
                    Output files: sample.ca.R1.fastq, sample.ca.R2.fastq, sample.ca.log. 
                    """
                    )
    parser.add_argument("sample_name", help="sample name - used as output prefix")

    parser.add_argument("fastq_F_file", help="fastq file with forward reads (.fastq or .fastq.gz)")

    parser.add_argument("fastq_R_file", nargs='?', default=None, help="Optional: fastq file with reverse reads (.fastq or .fastq.gz)")

    parser.add_argument("-d", "--dry_run", action="store_true", default=False, help="print out commands and exit")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--trim_illumina", action="store_true", default=True, help="trim illumina adapters")

    parser.add_argument("--ampl", default="dop", help="Amplification protocol - used to remove specific primers. Possible values: dop, wga, none")

    parser.add_argument("--params", default="--trim-n --minimum-length 20", help="Additional parameters for cutadapt. \
        Default - trim terminal Ns and discard read pairs with at least one read shorter than 20 (bowtie2 seed length). \
        For WGA it is sometimes useful to increase error toleance with -e")

    return parser.parse_args()

def main(args):

    assert os.path.isfile(args.fastq_F_file)
    if args.fastq_R_file:
        assert os.path.isfile(args.fastq_R_file)
    
    if args.fastq_R_file == None:
        reads = [args.fastq_F_file]
    else:
        reads = [args.fastq_F_file, args.fastq_R_file]
    '''# plain or gzipped input
    if all(read.endswith('.fastq') for read in reads): 
        gzipped = False
    elif all(read.endswith('.fastq.gz') for read in reads):
        gzipped = True
    else:
        raise Exception('Improper read naming:\n%s\nValid file extensions are .fastq for \
         uncompressed and .fastq.gz for gzipped reads' % (args.fastq_F_file, args.fastq_R_file))
    '''
    ca_reads = [args.sample_name + '.ca.R1.fastq', args.sample_name + '.ca.R2.fastq']
    log_file = args.sample_name + '.ca.log'
    if not os.path.isfile(ca_reads[0]):
        # custom parameters
        cutadapt_opts = args.path_to_cutadapt + ' ' + args.params.strip("\'\"")
        # remove Illumna TruSeq adapters
        if args.trim_illumina:
            cutadapt_opts += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC '
            if len(reads) == 2:
               cutadapt_opts += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT '
        # remove primers from both ends of reads depending on the protocol
        if args.ampl == 'dop':
            cutadapt_opts += ' -n 3 -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG '
            if len(reads) == 2:
                cutadapt_opts += ' -A CCACATNNNNNNCTCGAGTCGG -G CCGACTCGAGNNNNNNATGTGG '
        elif args.ampl == 'wga':
            cutadapt_opts += '-n 3 -a CCAAACACACCCAACACAA -g TTGTGTTGGGTGTGTTTGG '
            if len(reads) == 2:
                cutadapt_opts += '-A CCAAACACACCCAACACAA -G TTGTGTTGGGTGTGTTTGG'
        elif args.ampl == 'none':
            cutadapt_opts += ' '
        else:
             raise Exception('Unknown amplification protocol. Known ones - dop, wga, none')
        # inputs and outputs. Note that reads left unpaired after trimming are discarded.
        if len(reads) == 1:
            cutadapt_opts += ' -o %s %s' % (ca_reads[0], reads[0])
        elif len(reads) == 2:
            cutadapt_opts += ' -o %s -p %s %s %s' % (ca_reads[0], ca_reads[1], reads[0], reads[1])

        sys.stderr.write(cutadapt_opts + ' 2> %s\n' % (log_file)) 
        if not args.dry_run:
            with open(log_file, 'w') as log:
                process = subprocess.Popen(cutadapt_opts.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
                (out, err) = process.communicate()
                log.write(out)
                sys.stderr.write(err)
                if process.returncode != 0:
                    sys.exit()
    else:
        sys.stderr.write('%s reads with removed adapters exist. OK!\n' % (args.sample_name + '.ca.R*.fastq'))

if __name__ == '__main__':
    main(parse_command_line_arguments())

