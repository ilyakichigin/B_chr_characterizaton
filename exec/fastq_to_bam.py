#!/usr/bin/env python

import subprocess
import sys
import os
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Aligns paired fastq files (preliminary processed with fastq_clean.py) to reference genome, prints commands used to stderr. 
                    Required program: bowtie2 (tested on v.2.1.0, 2.2.4). 
                    """
                    )

    parser.add_argument("sample", help="sample name - prefix of input files processed with cutadapt")

    parser.add_argument("-a", "--aligner", default="bt2", help="Aligner to use. Possible values: bt2 (bowtie2), bbm (bbmap)")

    parser.add_argument("-r", "--reference_genome", help="path to reference genome fasta. For bowtie, index files are expected in the same folder and with proper basename: e.g. canFam3 for canFam3.masked.fasta. For bbmap, index files are generated in the working directory.")

    parser.add_argument("-p", "--path_to_aligner", default="bowtie2", help="path to aligner program")

    parser.add_argument("-b", "--aligner_args", default="-p 1", help="Additional parameters for aligner specified as quoted string. For reference, see manuals.")


    return parser.parse_args()

def run_aligner(command, log_file):
    
    # run command logging to a file
    sys.stderr.write(command + ' 2> %s\n' % (log_file)) 
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    (out, err) = process.communicate()
    with open(log_file, 'w') as log:
        log.write(err)
    if process.returncode != 0:
        print "Something went wrong. Check %s for details." % (log_file)
        sys.exit()

def main(args):    

    if os.path.isfile(args.sample + '.ca.fastq'):
        single_end = True
        ca_fq = args.sample + '.ca.fastq'
    elif os.path.isfile(args.sample + '.ca.R1.fastq') and os.path.isfile(args.sample + '.ca.R2.fastq'):
        single_end = False
        f_ca_fq = args.sample + '.ca.R1.fastq'
        r_ca_fq = args.sample + '.ca.R2.fastq'
    
    reference = args.reference_genome.split('/')[-1]
    reference = reference.split('.')[0]
    
    out_sam = '%s.%s_%s.sam' % (args.sample,reference,args.aligner) # simplified name: 'ca' not included
    # alignment to reference - if not already done
    if (not os.path.isfile(out_sam)) or (os.path.getsize(out_sam) == 0):

        log_file = out_sam[:-3]+'log'

        if args.aligner == 'bt2':
            # assert that indexed genome exists
            reference_path = '/'.join(args.reference_genome.split('/')[:-1])+'/'+reference
            assert os.path.isfile(reference_path + '.1.bt2') # only first of six bowtie2 index files is checked
            if single_end:
                command = '%s %s -U %s -x %s -S %s' % (args.path_to_aligner, 
                    args.aligner_args, ca_fq, reference_path, out_sam)
            else:
                command = '%s %s -1 %s -2 %s -x %s -S %s' % (args.path_to_aligner, 
                    args.aligner_args, f_ca_fq, r_ca_fq, reference_path, out_sam)
        elif args.aligner == 'bbm':
            assert os.path.isfile(args.reference_genome)
            if single_end:
                command = '%s %s in%s ref=%s out=%s' % (args.path_to_aligner, 
                    args.aligner_args, ca_fq, args.reference_genome, out_sam)
            else:
                command = '%s %s in1=%s in2=%s ref=%s out=%s' % (args.path_to_aligner, 
                    args.aligner_args, f_ca_fq, r_ca_fq, args.reference_genome, out_sam)
        else:
            raise Exception('Invalid aligner!')
        run_aligner(command, log_file)

    else:
        print 'Alignment to reference genome exists. OK!'

if __name__ == '__main__':
    main(parse_command_line_arguments())
