#!/usr/bin/env python

import subprocess
import sys
import os
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes primers and Illumina adapter in PE mode, prints out programs used to stdout. 
                    Required program: cutadapt (tested on v.1.8)
                    1) Rename reads to include '/1' for forward and '/2' for reverse read (for cutadapt).
                    2) Trim Illumina TruSeq adapters, as well as DOP or WGA primers. Write log.
                    Output files: sample.ca.R1.fastq, sample.ca.R2.fastq, sample.ca.log. 
                    """
                    )
    parser.add_argument("fastq_F_file", help="fastq file with forward reads (.fastq)")

    parser.add_argument("fastq_R_file", help="fastq file with reverse reads (.fastq)")

    parser.add_argument("--sample_name", help="sample name - used as output prefix")

    parser.add_argument("-r", "--rename_only", action="store_true", help="only fix fastq format")

    parser.add_argument("-d", "--delimiter", default=" ", help="delimiter between cluster name and read number")

    parser.add_argument("--path_to_cutadapt", default="cutadapt", help="path to cutadapt binary")

    parser.add_argument("--ampl", default="dop", help="Amplification protocol - used to remove specific primers. Possible values: dop, wga, none")

    parser.add_argument("--params", default="--trim-n --minimum-length 20", help="Additional parameters for cutadapt. \
        Default - trim terminal Ns and discard read pairs with at least one read shorter than 20 (bowtie2 seed length). \
        For WGA it is sometimes useful to increase error toleance with -e")

    return parser.parse_args()
   
def rename_reads(in_file_name, delimiter):

    # This is done for cutadapt compatibility
    # Rename reads in fastq files - remove space between cluster name and number: 1 for F, 2 for R.
    #delimiter = delimiter.strip('\'\"')
    assert len(delimiter) == 1
    
    with open(in_file_name, 'rU') as in_file:
        in_file_base = in_file_name.split('/')[-1] # output to current folder
        out_file_name = '.'.join(in_file_base.split('.')[:-1]+['rn','fastq'])
        sys.stderr.write("Renaming reads in %s. Writing to %s.\n" % (in_file_name,out_file_name))
        with open(out_file_name, 'w') as out_file:
            for line in in_file:
                line_list = line.split(delimiter)
                if line.startswith('@') and len(line_list) > 1:
                    out_line = line_list[0] + '/' + line_list[1][0] + '\n' # leave \1 for F and \2 for R
                elif line.startswith('+'):
                    out_line = '+\n'
                else:
                    out_line = line
                out_file.write(out_line)

    return out_file_name

def main(args):
    
    assert args.fastq_F_file.endswith('.fastq') and args.fastq_R_file.endswith('.fastq') # do not accept gzipped and improperly named files
    # rename reads
    f_rn_fq = rename_reads(args.fastq_F_file, args.delimiter) # Filename returned. rn = renamed
    r_rn_fq = rename_reads(args.fastq_R_file, args.delimiter)
    # assign output filenames

    if not args.rename_only: 
        f_ca_fq = args.sample_name + '.ca.R1.fastq'
        r_ca_fq = args.sample_name + '.ca.R2.fastq'
        log_file = args.sample_name + '.ca.log'
        # rename and cutadapt if not already done
        if not os.path.isfile(f_ca_fq) and not os.path.isfile(r_ca_fq):
            
            # common options: remove Illumna TruSeq adapters, apply additional parameters.
            cutadapt_opts = args.path_to_cutadapt + ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                                                    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ' + \
                                                    args.params.strip("\'\"")
            # remove primers from both ends of reads depending on the protocol
            if args.ampl == 'dop':
                cutadapt_opts += ' -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG \
                                   -A CCACATNNNNNNCTCGAGTCGG -G CCGACTCGAGNNNNNNATGTGG -n 3'
            elif args.ampl == 'wga':
                cutadapt_opts += ' -a CCAAACACACCCAACACAA -g TTGTGTTGGGTGTGTTTGG \
                                   -A CCAAACACACCCAACACAA -G TTGTGTTGGGTGTGTTTGG -n 3'
            elif args.ampl == 'none':
                cutadapt_opts += ' '
            else:
                 raise Exception('Unknown amplification protocol. Known ones - dop, wga, none')
            # inputs and outputs. Note that reads left unpaired after trimming are discarded.
            cutadapt_opts += ' -o %s -p %s %s %s' % (f_ca_fq,r_ca_fq,f_rn_fq,r_rn_fq)
            with open(log_file, 'w') as log:
                sys.stderr.write(cutadapt_opts+'\n') 
                process = subprocess.Popen(cutadapt_opts.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
                (out, err) = process.communicate()
                log.write(out)
                sys.stderr.write(err)
                if process.returncode != 0:
                    sys.exit()
            sys.stderr.write('rm %s %s\n' % (f_rn_fq,r_rn_fq))
            os.remove(f_rn_fq)
            os.remove(r_rn_fq)
        else:
            print 'Files with removed adapters exist. OK!'
    else:
        print 'No adapter trimming. OK!'


if __name__ == '__main__':
    main(parse_command_line_arguments())

