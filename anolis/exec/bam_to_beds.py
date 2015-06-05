#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    From bam file create reads and positions bed files. 
                    Requires: bedtools (tested on v.2.17.0)
                    """
                    )
    parser.add_argument("bam_file", help="input bam file")
    parser.add_argument("--path_to_bedtools", default="bedtools",help="path to bedtools")

    return parser.parse_args()

def run_bedtools(bam_file, path_to_bedtools):

    # Using bedtools generate .bed files with reads and positions. Returns output file names.
    
    read_bed_name = bam_file[:-4]+'.reads.bed'
    pos_bed_name = bam_file[:-4]+'.pos.bed'   
    
    command_list = [ # [command string, output file name]
            [path_to_bedtools + ' bamtobed -i ' + bam_file, read_bed_name],
            [path_to_bedtools + ' merge -n -i ' + read_bed_name, pos_bed_name]
            ]    

    for command in command_list:
        sys.stderr.write(command[0]+'\n')
        process = subprocess.Popen(command[0].split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        process.wait()
        if err:
            sys.stderr.write(err)
            sys.exit()
        with open(command[1],'w') as out_file:
            sys.stderr.write('Writing output to file '+command[1]+'\n')
            out_file.write(out)
        
        
if __name__ == '__main__':
    args = parse_command_line_arguments()
    assert args.bam_file.endswith('bam')
    wg_bed_files = run_bedtools(args.bam_file, args.path_to_bedtools)
