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
    parser.add_argument("--bedtools_path", default="bedtools",help="path to bedtools")
    parser.add_argument("-d", "--dry_run", action="store_true", 
                        default=False, help="print out commands and exit")
    return parser.parse_args()

def run_command(command, verbose=True, dry_run=False, use_shell=False):
    if verbose:
        sys.stderr.write(command)
    if not dry_run:
        if not use_shell:
            command = command.split(' ')
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=use_shell)
        (out, err) = p.communicate()
        if err:
            sys.stderr.write(err)
            sys.exit()

        return out

def main(args):
    if args.bam_file.endswith('.filter.bam'):
        pos_bed_file = args.bam_file[:-11] + '.pos.bed'
    elif args.bam_file.endswith('.bam'):
        pos_bed_file = args.bam_file[:-4]+'.pos.bed'
    else:
        raise Exception('Invalid BAM file suffix:%s\nAcceptable values: *.filter.bam *.bam\n')
    
    # version-dependent merge
    if not os.path.isfile(pos_bed_file):
        bedtools_version_string = run_command(args.bedtools_path + ' --version', verbose = False)
        bedtools_version = bedtools_version_string.split()[1].split('.')[1]
        if int(bedtools_version) >= 20:
            bedtools_merge_command = args.bedtools_path + ' merge -c 1 -o count'
        else:
            bedtools_merge_command = args.bedtools_path + ' merge -n -i'

        command = '%s bamtobed -i %s | %s sort | %s ' % (args.bedtools_path, args.bam_file, 
            args.bedtools_path, bedtools_merge_command)
        pos_bed = run_command(command, verbose=True, dry_run=args.dry_run, use_shell=True)
        sys.stderr.write('> %s\n' % pos_bed_file)
        if not args.dry_run:
            with open(pos_bed_file, 'w') as outf:
                outf.write(pos_bed)
    else:
        sys.stderr.write('%s BED file with read position exists. OK!\n' % pos_bed_file)

if __name__ == '__main__':
    main(parse_command_line_arguments())
    

