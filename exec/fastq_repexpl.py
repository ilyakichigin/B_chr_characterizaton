#!/usr/bin/env python

import sys
import argparse

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Converts paired reads from Illumina 1.8+ fastq processed with cutadapt to interleaved fasta. 
                    Suitable as input for RepeatExplorer with 'All sequence reads are paired' option.
                    """
                    )
    #parser.add_argument("fastq_F_file", help="fastq file with forward reads (.fastq)")
    #parser.add_argument("fastq_R_file", help="fastq file with reverse reads (.fastq)")
    parser.add_argument("cutadapt_prefix", help="prefix of fastq files processed with cutadapt. Format <prefix>.<F|R.ca.fastq>")    
    parser.add_argument("-r","--rename", action="store_true",help="rename reads to numeric")
    return parser.parse_args()
   
def fastq_to_re_fasta(args):

    # From two fastq files generate interleaved fasta file with added '/1' for forward and '/2' for reverse read
    
    f_fastq = args.cutadapt_prefix + '.ca.R1.fastq'
    r_fastq = args.cutadapt_prefix + '.ca.R2.fastq'
    o_fasta = args.cutadapt_prefix + '.ca.re.fasta'
    rename = args.rename
    
    min_read_length = 18 # minimum read length defined as word length in megablast at cluster annotation
    i=0 # read number
    k=0 # line number
    fname = ''
    rname = ''
    print 'Output file: ' + o_fasta
    with open(f_fastq, 'rU') as f_file, open(r_fastq, 'rU') as r_file, open(o_fasta, 'w') as out:
        for fline in f_file:
            rline = r_file.next()
            if fname != '' and rname != '': # read sequence - next line after name
                if len(fline) > min_read_length + 1 and len(rline) > min_read_length + 1: # exclude pairs with short reads
                    out.write(fname+fline+rname+rline)
                fname = ''
                rname = ''
            elif fline.startswith('@') and k%4 == 0: # get read name
                i+=1
                if rename:
                    fname='>%d/1\n'%i
                    rname='>%d/2\n'%i
                else:
                    fname='>'+fline[1:-1]+'/1\n'
                    rname='>'+rline[1:-1]+'/2\n'
            k+=1

if __name__ == '__main__':
    fastq_to_re_fasta(parse_command_line_arguments())

    
