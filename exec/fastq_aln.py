#!/usr/bin/env python

import subprocess
import sys
import os
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Aligns fastq files (preliminary processed with fastq_clean.py) to reference genome 
                    Required program: bowtie2 (tested on v.2.1.0, 2.2.4). 
                    """
                    )

    parser.add_argument("fastq_F_file", help="fastq file with forward reads")

    parser.add_argument("fastq_R_file", nargs='?', default=None, help="Optional: fastq file with reverse reads")

    parser.add_argument("-a", "--aligner", default="bt2", 
        help="aligner to use. Possible values: bt2 (bowtie2), bwa, bbm (bbmap)")

    parser.add_argument("-r", "--reference_genome", help="path to reference genome fasta")

    parser.add_argument("-p", "--path_to_aligner", default="bowtie2", help="path to aligner program")

    parser.add_argument("-b", "--aligner_args", default="-p 1", 
        help="additional parameters for aligner specified as quoted string. For reference, see manuals.")

    parser.add_argument("-d", "--dry_run", action="store_true", default=False, 
        help="print out commands and exit")


    return parser.parse_args()

def run_aligner(command, log_file, dry_run=False, stdout_file=None):
    
    # run command logging to a file
    sys.stderr.write(command)
    if stdout_file is not None:
        sys.stderr.write(' > %s' % (stdout_file))
    sys.stderr.write(' 2> %s\n' % (log_file)) 
    if not dry_run:
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        (out, err) = process.communicate()
        with open(log_file, 'w') as log:
            log.write(err)
        if process.returncode != 0:
            raise Exception("Something went wrong. Check %s for details." % (log_file))
        if stdout_file is not None:
            with open(stdout_file, 'w') as out_sam:
                out_sam.write(out)

def main(args):    

    # outputs
    sample = args.fastq_F_file.split('.')[0]
    rg = args.reference_genome.split('/')[-1].split('.')[0]
    rg_path = '/'.join(args.reference_genome.split('/')[:-1]) + '/' + rg
    out_sam = '%s.%s.sam' % (sample, rg) 
    log_file = '%s.%s.log' % (sample, rg) 
    
    if (not os.path.isfile(out_sam)) or (os.path.getsize(out_sam) == 0):

        if args.aligner == 'bt2':
            # bowtie2 prefix example: canFam3 will be used as rg for canFam3.masked.fasta
            if not os.path.isfile(rg_path + '.1.bt2'):
                raise Exception('bowtie2 index not found. Expected path and prefix: %s' % rg_path)
            if args.fastq_R_file == None:
                command = '%s %s -U %s -x %s -S %s' % (args.path_to_aligner, 
                    args.aligner_args, args.fastq_F_file, rg_path, out_sam)
            else:
                command = '%s %s -1 %s -2 %s -x %s -S %s' % (args.path_to_aligner, 
                    args.aligner_args, args.fastq_F_file, args.fastq_R_file, rg_path, out_sam)
            run_aligner(command, log_file, args.dry_run)

        elif args.aligner == 'bwa':
            if args.path_to_aligner == "bowtie2":
                args.path_to_aligner = "bwa"
            if not os.path.isfile(rg_path + '.bwt'):
                raise Exception('bowtie2 index not found. Expected path and prefix: %s' % rg_path)
            if args.fastq_R_file == None:
                command = '%s %s %s %s' % (args.path_to_aligner, 
                    args.aligner_args, rg_path, args.fastq_F_file)
            else:
                command = '%s %s %s %s %s' % (args.path_to_aligner, 
                    args.aligner_args, rg_path, args.fastq_F_file, args.fastq_R_file)
            run_aligner(command, log_file, args.dry_run, out_sam)

        elif args.aligner == 'bbm':
            if args.path_to_aligner == "bowtie2":
                args.path_to_aligner = "bbmap.sh"
            assert os.path.isfile(args.reference_genome)
            if args.fastq_R_file == None:
                command = '%s %s in%s ref=%s out=%s' % (args.path_to_aligner, 
                    args.aligner_args, args.fastq_F_file, args.reference_genome, out_sam)
            else:
                command = '%s %s in1=%s in2=%s ref=%s out=%s' % (args.path_to_aligner, 
                    args.aligner_args, args.fastq_F_file, args.fastq_R_file, args.reference_genome, out_sam)
            run_aligner(command, log_file, args.dry_run)

        else:
            raise Exception('Invalid aligner!')
        

    else:
        sys.stderr.write('%s alignment to reference genome exists. OK!\n' % out_sam)

if __name__ == '__main__':
    main(parse_command_line_arguments())
