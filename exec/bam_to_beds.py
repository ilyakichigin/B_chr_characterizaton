#!/usr/bin/env python


import subprocess
import sys
import argparse
import os


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
    srt_read_bed_name = bam_file[:-4]+'.srt.reads.bed'
    pos_bed_name = bam_file[:-4]+'.pos.bed'   
    
    command_list = [ # [command string, output file name]
            [path_to_bedtools + ' bamtobed -i ' + bam_file, read_bed_name],
            [path_to_bedtools + ' sort -i ' + read_bed_name, srt_read_bed_name],
            [path_to_bedtools + ' merge -n -i ' + srt_read_bed_name, pos_bed_name]
            ]    

    for command in command_list:
        sys.stderr.write('%s > %s\n' % (command[0],command[1]))
        process = subprocess.Popen(command[0].split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if err:
            sys.stderr.write(err)
            sys.exit()
        with open(command[1],'w') as out_file:
            out_file.write(out)
    sys.stderr.write('mv %s %s\n' % (srt_read_bed_name,read_bed_name))    
    os.rename(srt_read_bed_name,read_bed_name)

def main(args):
    assert args.bam_file.endswith('bam')
    wg_bed_files = run_bedtools(args.bam_file, args.path_to_bedtools)
        
if __name__ == '__main__':
    main(parse_command_line_arguments())
    

