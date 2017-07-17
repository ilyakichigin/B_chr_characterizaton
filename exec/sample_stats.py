#!/usr/bin/env python

import subprocess
import tempfile
import argparse
import pysam
import sys
import os

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """Calculates statistics of sample sequencing and alignment to reference genome.
                    """
                    )
    parser.add_argument("sample", help="sample name - used as prefix")
    parser.add_argument("genome", help="target genome prefix")
    parser.add_argument("-f","--f_reads", default='', help="forward reads .fastq[.gz]")
    parser.add_argument("-r","--r_reads", default='', help="reverse reads .fastq[.gz]")
    parser.add_argument("-t", "--target_regions", default=None, help="verified target regions .bed")
    parser.add_argument("-o", "--output", default=None, help="Output filename. default <sample>.stats.txt")

    return parser.parse_args()

def generate_filenames(sample, genome):

    names = dict()
    names['ca_log'] = '.'.join([sample,'ca','log'])
    names['filter_log'] = '.'.join([sample,genome,'filter','log'])
    names['pos_bed'] = '.'.join([sample,genome,'pos','bed'])

    return names

def bed_stats(f):

    i = 0
    pos_size = 0
    read_num = 0

    for line in f:
        ll = line.split('\t')
        i += 1
        pos_size += int(ll[2])-int(ll[1])
        if len(ll) > 3:
            read_num += int(ll[3])

    return(str(i), str(pos_size), str(float(read_num)/i), str(float(pos_size)/i))

def del_1000_sep(s):
    return ''.join(s.split(','))

def main(args):
    
    names = generate_filenames(args.sample, args.genome)
    # collect stats as list of tuples
    stats = [('sample',args.sample)]

    # fastq stats
    with open(names['ca_log']) as f:
        for line in f:
            ll = line.split()
            if line.startswith('Total read pairs processed'):
                stats.append(('Total read pairs', del_1000_sep(ll[-1])))
            elif line.startswith('Pairs written (passing filters)'):
                stats.append(('Trimmed read pairs',del_1000_sep(ll[-2])))
            elif line.startswith('Total reads processed'):
                stats.append(('Total reads',del_1000_sep(ll[-1])))
            elif line.startswith('Reads written (passing filters)'):
                stats.append(('Trimmed reads',del_1000_sep(ll[-2])))
            elif line.startswith('Total basepairs processed'):
                stats.append(('Total bp',del_1000_sep(ll[-2])))
            elif line.startswith('Total written (filtered):'):
                stats.append(('Trimmed bp',del_1000_sep(ll[-3])))

    with open(names['filter_log']) as f:
        for line in f:
            if line.startswith('Target aligned fragments'):
                stats.append(('Aligned reads', line.split()[-1]))
            elif line.startswith('Filtered Q'):
                stats.append(('Filtered reads', line.split()[-1]))
            elif line.startswith('Contamination Q'):
                stats.append(('Contamination reads', line.split()[-1]))

    # pos file stats
    with open(names['pos_bed']) as f:
        pos_stat = bed_stats(f)
        stats += ('Read positions', pos_stat[0]),('Position sum bp', pos_stat[1])
        stats += ('Reads per position', pos_stat[2]),('Position mean bp', pos_stat[3])
    # regpos stats
    if args.target_regions is not None:
        with open(args.target_regions) as reg_bed:
            reg_stat = bed_stats(reg_bed)
            stats += ('Regions', reg_stat[0]),('Regions bp', reg_stat[1])

        with tempfile.TemporaryFile() as regpos_bed:
            bt_command = ['bedtools','intersect','-a',args.target_regions,'-b',names['pos_bed']]
            p = subprocess.Popen(bt_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = p.communicate()
            if p.returncode != 0:
                raise Exception("bedtools intersect failed:\n\n%s" % (err))
            regpos_bed.write(out)
            regpos_bed.seek(0)
            regpos_stat = bed_stats(regpos_bed)
            stats += ('Read positions in regions', regpos_stat[0]),('Position sum bp in regions', regpos_stat[1])
            stats += ('Position in region/region bp', str(float(regpos_stat[1])/int(reg_stat[1])))
            stats += ('Position in region/genome positions bp', str(float(regpos_stat[1])/int(pos_stat[1])))
    #output
    if args.output is None:
        args.output = args.sample + '.stats.txt'
    with open(args.output,'w') as f:
        f.write('\n'.join(['\t'.join(s) for s in stats]))

if __name__ == '__main__':
    args = parse_command_line_arguments()
    main(args)
