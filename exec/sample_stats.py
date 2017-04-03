#!/usr/bin/env python

import subprocess
import tempfile
import argparse
import pysam
import gzip
import sys
import re
import os

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """Calculates statistics of sample sequencing and alignment to reference genome.
                    """
                    )
    parser.add_argument("sample", help="sample name - used as prefix")
    parser.add_argument("F_reads", help="forward reads .fastq")
    parser.add_argument("R_reads", help="reverse reads .fastq")
    parser.add_argument("-r", "--regions", default=None, help="verified target regions .bed")

    return parser.parse_args()

def get_filenames(sample):

    names = dict()
    names['f_trim_fq'] = sample + '.ca.R1.fastq'
    names['r_trim_fq'] = sample + '.ca.R2.fastq'

    fb_names = [f for f in os.listdir('.') if re.match(r'.*\.filter\.bam$', f)]
    if len(fb_names) < 2:
        names['filtered_bam'] = ''.join(fb_names)
    else:
        raise Exception('Too many *filter.bam files in current folder %s' % str(fb_names))

    cb_names = [f for f in os.listdir('.') if re.match(r'.*\.contam\.bam$', f)]
    if len(cb_names) < 2:
        names['contam_bam'] = ''.join(cb_names)
    else:
        raise Exception('Too many *contam.bam files in current folder %s' % str(cb_names))

    pb_names = [f for f in os.listdir('.') if re.match(r'.*\.pos\.bed$', f)]
    if len(pb_names) < 2:
        names['pos_bed'] = ''.join(pb_names[0])
    else:
        raise Exception('Too many *pos.bed files in current folder')

    return names

def fastq_stats(fastq_name):

    i = 0
    readlen = 0
    
    if fastq_name.endswith('.fastq'):
        f = open(fastq_name)
    elif fastq_name.endswith('.fastq.gz'):
        f = gzip.open(fastq_name)
    else:
        raise Exception('Improper read naming: %s' % (fastq_name))

    for line in f:
        i += 1
        if i%4 == 1 and not line.startswith('@'):
            raise Exception('Corrupt read file: %s' % (fastq_name))
        if i%4 == 2:
            readlen += len(line.strip('\n'))
    f.close()

    return [i/4, readlen]

def bam_stats(filename):

    n_reads = 0
    if filename != '':
        # from https://www.biostars.org/p/1890/ - get number of mapped reads
        for l in pysam.idxstats(filename).split('\n'):
            ll = l.rstrip('\n').split('\t')
            if len(ll) == 4:
                n_reads += int(ll[2])
        return n_reads
    else:
        return 'NA'

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

    return(i, pos_size, float(read_num)/i)

def main(args):
    
    names = get_filenames(args.sample)
    stats = [('sample',args.sample)]

    # init fastq stats
    f_init_stat = fastq_stats(args.F_reads)
    r_init_stat = fastq_stats(args.F_reads)
    stats += ('Initial reads', str(f_init_stat[0] + r_init_stat[0])), ('Initial bp', str(f_init_stat[1] + r_init_stat[1]))
    # trimmed fastq stats
    f_trim_stat = fastq_stats(names['f_trim_fq'])
    r_trim_stat = fastq_stats(names['r_trim_fq'])
    stats += ('Trimmed reads', str(f_trim_stat[0] + r_trim_stat[0])), ('Trimmed bp', str(f_trim_stat[1] + r_trim_stat[1]))

    # post-filtering bam stats
    stats.append(('Contamination reads', str(bam_stats(names['contam_bam']))))
    stats.append(('Target reads after filtering', str(bam_stats(names['filtered_bam']))))
    # pos file stats
    with open(names['pos_bed']) as f:
        pos_stat = bed_stats(f)
        stats += ('Read positions', str(pos_stat[0])),('Position sum bp', str(pos_stat[1])),('Reads per position', str(pos_stat[2]))
    stats.append(('Position mean bp', str(float(pos_stat[1])/pos_stat[0])))
    # regpos stats
    if args.regions is not None:
        with open(args.regions) as reg_bed:
            reg_stat = bed_stats(reg_bed)
            stats += ('Regions', str(reg_stat[0])),('Regions bp', str(reg_stat[1]))

        with tempfile.TemporaryFile() as regpos_bed:
            bt_command = ['bedtools','intersect','-a',args.regions,'-b',names['pos_bed']]
            p = subprocess.Popen(bt_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = p.communicate()
            if p.returncode != 0:
                raise Exception("bedtools intersect failed:\n\n%s" % (err))
            regpos_bed.write(out)
            regpos_bed.seek(0)
            regpos_stat = bed_stats(regpos_bed)
            stats += ('Read positions in regions', str(regpos_stat[0])),('Position sum bp in regions', str(regpos_stat[1]))
        stats.append(('Position in region/region bp', str(float(regpos_stat[1])/reg_stat[1])))
        stats.append(('Position in region/genome positions bp', str(float(regpos_stat[1])/pos_stat[1])))
    print '\n'.join(['\t'.join(s) for s in stats])

if __name__ == '__main__':
    args = parse_command_line_arguments()
    main(args)
