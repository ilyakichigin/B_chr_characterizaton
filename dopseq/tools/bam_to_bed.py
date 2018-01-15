#!/usr/bin/env python


import os
import sys
import argparse

from dopseq.tools import utils


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Convert BAM to BED, merge overlapping reads into 
                    read positions
                    """
                    )
    parser.add_argument("in_bam", help="input bam file")
    parser.add_argument("--out_bed",
                        default="test.bed",
                        help="output BED file")
    parser.add_argument("--bedtools_path",
                        default="bedtools",
                        help="path to bedtools")
    parser.add_argument("-d", "--dry_run", 
                        action="store_true", 
                        default=False, 
                        help="print out commands and exit")
    return parser.parse_args()


def bedtools_release(bedtools_path, dry_run):
    """Given path to bedtools return its release number 
    (second number in version string)
    """
    bedtools_version_command = bedtools_path + ' --version'
    bedtools_version = utils.run_command(bedtools_version_command, 
                                         verbose=False,
                                         dry_run=dry_run,
                                         return_out=True)
    bedtools_release = bedtools_version.split()[1].split('.')[1]

    return int(bedtools_release)


def main(args):
    
    init_bed = utils.run_command('%s bamtobed -i %s' % (args.bedtools_path, 
                                                        args.in_bam), 
                                 verbose=True, dry_run=args.dry_run, 
                                 return_out=True)

    sort_bed = utils.run_command('%s sort -i -' % args.bedtools_path, 
                                 verbose=True, dry_run=args.dry_run, 
                                 stdin = init_bed, return_out=True)

    utils.run_command('%s merge -c 1 -o count -i -' % args.bedtools_path, 
                      verbose=True, dry_run=args.dry_run, 
                      stdin = sort_bed, outfile=args.out_bed)

    
if __name__ == '__main__':
    main(parse_command_line_arguments())
    

