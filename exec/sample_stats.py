#!/usr/bin/env python

import argparse
import pysam

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """Calculates statistics of args.sample sequence and alignment to reference genome.
                    """
                    )
    parser.add_argument("args.sample", nargs="*", help="args.sample name - used as prefix for ")
    parser.add_argument("args.F_reads", nargs="*", help="file with forward reads")
    parser.add_argument("args.R_reads", nargs="*", help="file with reverse reads")
    parser.add_argument("args.t_genome", nargs="*", help="target genome prefix")
    parser.add_argument("args.c_genome", nargs="*", help="contamination genome prefix")
    parser.add_argument("-H", "--header", required=False, action="store_true",
                        help="print out header and exit")
    #parser.add_argument("-tr", "--total_reads", type=int, default=70,
    #                help="minimum number of reads in whole region needed to consider it a true region")

    return parser.parse_args()

def fastq_stats(fastq_name):
    i = 0
    readlen = 0
    with open(fastq_name) as f:
        for line in f:
            i += 1
            if i%4 == 2:
                readlen += len(line.strip('\n'))

        return [i/4, readlen]

def bam_stats(filename):
    # from https://www.biostars.org/p/1890/ - number of mapped reads
    return reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(filename) ])

def main(args):
    
    stats = [args.sample, 0, 0, 0, 0]
    # init fastq stats
    for fname in (args.F_reads, args.R_reads):
        fstat = fastq_stats(fname)
        stats[1] += fstat[0] # reads
        stats[2] += fstat[1] # bp
    # cutadapt fastq stats
    for direction in ('F','R'):
        fname = '.'.join([args.sample,direction,'ca','fastq'])
        fstat = fastq_stats(fname)
        stats[3] += fstat[0] # reads
        stats[4] += fstat[1] # bp
    # number of reads mapped to contamination genome
    stats.append(bam_stats('.'.join([args.sample,args.c_genome,'contam.bam'])))
    stats.append(bam_stats('.'.join([args.sample,args.t_genome,'filter.bam'])))
    # pos file stats
    pos_file = '.'.join([args.sample,args.t_genome,'filter.pos.bed'])
    i = 0
    pos_size = 0
    l = []
    with open(pos_file) as f:
        for line in f:
            ll = line.split('\t')
            i += 1
            pos_size += int(ll[2])-int(ll[1])
            l.append(int(ll[3]))
    stats.extend([i,pos_size,sum(l)/float(len(l))])

    stats = [str(i) for i in stats]
    with open(args.sample + '.stats.txt', 'w') as outf:
        outf.write( '\t'.join(stats)+'\n' )

if __name__ == '__main__':
    args = parse_command_line_arguments()
    if args.header:
        print 'init_reads\tinit_bp\ttrimmed_reads\ttrimmed_bp\tcontam_reads\ttarget_reads\tpositions\tpos_bp\tmean_coverage'
    else:
        main(args)
