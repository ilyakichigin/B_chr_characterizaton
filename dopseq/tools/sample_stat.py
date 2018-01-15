#!/usr/bin/env python

import subprocess
import tempfile
import argparse
import pysam
import sys
import os
from collections import OrderedDict

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """Calculates statistics of sample sequencing and alignment to reference genome.
                    """
                    )
    parser.add_argument("ca_log", help="cutadapt log")
    parser.add_argument("filter_log", help="contam_filter log")
    parser.add_argument("pos_bed", 
                        help="BED file with merged mapped read positions")
    parser.add_argument("-o", "--output", default=None, 
                        help="Output filename. default <sample>.stats.txt")
    parser.add_argument("-d", "--dry_run", default=False, 
                        help="Perform only a dry run")

    return parser.parse_args()


def bed_stats(f):
    """Given open file in BED format return
    - number of intervals
    - total size of bed positions
    - mean number of reads (column 4) per position
    - mean position size
    """
    i = 0
    pos_size = 0
    read_num = 0

    for line in f:
        ll = line.split('\t')
        i += 1
        pos_size += int(ll[2])-int(ll[1])
        if len(ll) > 3:
            read_num += int(ll[3])

    reads_per_position = (0 if i == 0 else float(read_num)/i)
    mean_position_size = (0 if i == 0 else float(pos_size)/i)

    return (i, pos_size, reads_per_position, mean_position_size)


def del_1000_sep(s):
    """Remove thousand separator from string representation of number.
    """
    return int(''.join(s.split(',')))


def main(args):

    sample = args.pos_bed.split('.')[0]
    if args.output is None:
        args.output = sample + '.stats.txt'    

    sys.stderr.write('Output file: %s\n' % args.output)

    if not args.dry_run:
        # collect stats as list of tuples
        stats = OrderedDict()
        stats['sample'] = sample

        # fastq stats
        with open(args.ca_log) as f:
            for line in f:
                ll = line.split()
                if line.startswith('Total read pairs processed'):
                    stats['Total reads'] = del_1000_sep(ll[-1]) * 2
                elif line.startswith('Pairs written (passing filters)'):
                    stats['Trimmed reads'] = del_1000_sep(ll[-2]) * 2
                elif line.startswith('Total reads processed'):
                    stats['Total reads'] = del_1000_sep(ll[-1]) 
                elif line.startswith('Reads written (passing filters)'):
                    stats['Trimmed reads'] = del_1000_sep(ll[-2])
                elif line.startswith('Total basepairs processed'):
                    stats['Total bp'] = del_1000_sep(ll[-2])
                elif line.startswith('Total written (filtered):'):
                    stats['Trimmed bp'] = del_1000_sep(ll[-3])

        with open(args.filter_log) as f:
            for line in f:
                if line.endswith('segments in target alignment\n'):
                    stats['Reads aligned'] = line.split()[0]
                elif line.endswith('belonging to contaminant genome\n'):
                    stats['Reads contaminant'] = line.split()[0]
                elif line.endswith('discard unmapped)\n'):
                    stats['Reads not passing quality filters'] = line.split()[0]
                elif line.endswith('and filtering\n'):
                    stats['Reads retained'] = line.split()[0]
                
        # pos file stats
        with open(args.pos_bed) as f:
            pos_stat = bed_stats(f)
            stats['Merged read positions'] =  pos_stat[0]
            stats['Position total bp'] = pos_stat[1]
            stats['Position mean bp'] = '%.2f' % pos_stat[3]
            stats['Reads per position'] = '%.2f' % pos_stat[2]
            
        with open(args.output,'w') as out:
            for key in stats:
                out.write('%s\t%s\n' % (key, str(stats[key])))
        

if __name__ == '__main__':
    main(parse_command_line_arguments())
