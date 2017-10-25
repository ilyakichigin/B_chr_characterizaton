#!/usr/bin/env python

import os
import sys
import pysam
import tempfile
import argparse

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes contamination reads by comparing MAPQs to target and cotamination genomes.
                    Input: sam/bam alignments to target and contamination genomes.
                    Output: bam alignment: sample.genome.filter.bam
                    Pysam python package and samtools are required \
                    (tested on v.0.8.1 and v.0.1.19-44428cd, respectively)
                    """
                    )
    parser.add_argument("target_bam",
                        help="BAM alignment to target genome")
    parser.add_argument("contam_bam",
                        help="BAM alignment to contamination genome")
    parser.add_argument("-m", '--min_quality', default=20,
                        help="Minimum quality for filtered file. Default: 20.")
    parser.add_argument("-l", '--min_length', default=20,
                        help="Minimum alignment length for filtered file. Default: 20 bp.")
    parser.add_argument("-d", "--dry_run", action="store_true", 
                        default=False, help="print out commands and exit")

    return parser.parse_args()


def bam_sort(in_file, out_file, name_sort=False):
    samtools_rel = pysam.version.__samtools_version__[0] 
    if samtools_rel == '1':
        if name_sort:
            pysam.sort('-n', '-T', '/tmp/bam_nsort', '-o', out_file, in_file)
        else:
            pysam.sort('-T', '/tmp/bam_sort', '-o', out_file, in_file)
    elif samtools_rel == '0':
        raise Exception('Unsupported samtools verion: %s. Please update to pysam with samtools version 1.*.*.' % (pysam.version.__samtools_version__[0]))
    else:
        raise Exception('Unrecognized samtools verion: %s. Supported versions: 1.*.*.' % (pysam.version.__samtools_version__[0]))

def compare_mapqs(treads, creads, filter_file, min_qual, min_len, i, j):
    if sum([x.mapping_quality for x in treads])/len(treads) >= sum([x.mapping_quality for x in creads])/len(creads):
        for read in treads:
            if read.mapping_quality >= min_qual and read.reference_length >= min_len: # filter
                i += 1
                filter_file.write(read)
    else:
        for read in creads:
            if read.mapping_quality >= min_qual and read.reference_length >= min_len: # contam
                j += 1

    return (i, j)

def main(args):
    # output: filtered bam alignment
    filter_bam = args.target_bam[:-4] + '.filter.bam'
    log_file = args.target_bam[:-4] + '.filter.log'
    if not os.path.isfile(filter_bam):
        sys.stderr.write('%s -m %d %s %s > %s 2> %s\n' 
                % (os.path.realpath(__file__), args.min_quality, 
                args.target_bam, args.contam_bam, filter_bam, log_file))
        if not args.dry_run:
            min_qual = int(args.min_quality)
            min_len = int(args.min_length)
            
            # name-sort inputs
            t = tempfile.NamedTemporaryFile(suffix = '_t.bam', delete=False)
            c = tempfile.NamedTemporaryFile(suffix = '_c.bam', delete=False)
            tname = t.name
            cname = c.name
            t.close()
            c.close()
            bam_sort(args.target_bam, tname, name_sort=True)
            bam_sort(args.contam_bam, cname, name_sort=True)

            # filter to the intermediate output
            f = tempfile.NamedTemporaryFile(suffix = '_f.bam', delete=False)
            fname = f.name
            f.close()
            with pysam.AlignmentFile(tname) as tfile, pysam.AlignmentFile(cname) as cfile:
                filter_file = pysam.AlignmentFile(fname, 'wb', template=tfile)
                tf = 0; cf = 0; i = 0; j = 0
                tread = tfile.next()
                treads = [tread]
                cread = cfile.next()
                creads = [cread]
                if tread.query_name != cread.query_name:
                    raise Exception('First read names in name-sorted BAMs do not match')
                prev_query = tread.query_name
                for tread in tfile:
                    if not tread.is_unmapped:
                        tf += 1
                    if tread.query_name == prev_query: # accumulate targets
                        treads.append(tread)
                    else: # analyze contam
                        for cread in cfile:
                            if not cread.is_unmapped:
                                cf += 1
                            if cread.query_name == prev_query: # accumulate contam
                                    creads.append(cread)
                            else:
                                break
                        (i, j) = compare_mapqs(treads, creads, filter_file, min_qual, min_len, i, j)
                        treads = [tread]
                        creads = [cread]
                        prev_query = tread.query_name
                for cread in cfile:
                    cf += 1
                    if cread.query_name == prev_query: # accumulate contam
                        creads.append(cread)
                    else:
                        raise Exception('%s file has read %s not present in %s' % (args.contam_bam, cread.query_name, args.target_bam))
                (i, j) = compare_mapqs(treads, creads, filter_file, min_qual, min_len, i, j)
                filter_file.close()
                with open(log_file, 'w') as outfile:
                    outfile.write('Target aligned fragments:\t%d\n' % (tf))
                    outfile.write('Contam aligned fragments:\t%d\n' % (cf))
                    outfile.write('Contamination Q%d fragments (>=%d bp):\t%d\n' % (min_qual, min_len, j))
                    outfile.write('Filtered Q%d fragments (>=%d bp):\t%d\n' % (min_qual, min_len, i))
                #sys.stderr.write('%d init target fragments %d init contam fragments %d filtered %d contam\n' % (tf, cf, i, j))

            # coordinate-sort output
            bam_sort(fname, filter_bam)
            #sys.stderr.write('Writing %d alignments to filtered output file %s.\n'%(i, filter_bam))
            #sys.stdout.write('3 output files: filtered %s, contamination %s, unmapped %s\n'%(filter_bam, contam_bam, unmap_bam))

            #clean-up
            os.unlink(tname)
            os.unlink(cname)
            os.unlink(fname)
    else:
        sys.stderr.write('%s filtered alignment exists. OK!\n' % filter_bam)

if __name__ == '__main__':
    main(parse_command_line_arguments())
    
