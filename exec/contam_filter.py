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
    parser.add_argument("target_sam",
                        help="SAM alignment to target genome")
    parser.add_argument("contam_sam",
                        help="SAM alignment to contamination genome")
    parser.add_argument("-m", '--min_quality', default=20,
                        help="Minimum quality for filtered file. Default: 20.")
    parser.add_argument("-d", "--dry_run", action="store_true", 
                        default=False, help="print out commands and exit")

    return parser.parse_args()


def bam_sort(in_file, out_file, name_sort=False):
    samtools_rel = pysam.version.__samtools_version__[0] 
    if samtools_rel == '0':
        if name_sort:
            pysam.sort('-n', in_file, out_file[:-4])
        else:
            pysam.sort(in_file, out_file[:-4])
    elif samtools_rel == '1':
        if name_sort:
            pysam.sort('-n', '-T', '/tmp/bam_nsort', '-o', out_file, in_file)
        else:
            pysam.sort('-T', '/tmp/bam_sort', '-o', out_file, in_file)
    else:
        raise Exception('Invalid samtools verion: %s. Valid versions 0.*.*, 1.*.*.' % (pysam.version.__samtools_version__[0]))

def main(args):

    # output: filtered bam alignment
    filter_bam = args.target_sam[:-4] + '.filter.bam'
    if not os.path.isfile(filter_bam):
        if args.dry_run:
            sys.stderr.write('Using pysam version %s with samtools version %s.\nInput files: %s %s Output file: %s\n' 
                % (pysam.__version__, pysam.version.__samtools_version__, 
                    args.target_sam, args.contam_sam, filter_bam))
        else:
            sys.stderr.write('Target: %s. Contamination: %s.\n'%(args.target_sam, args.contam_sam))
            min_qual = int(args.min_quality)
            
            # name-sort inputs
            t = tempfile.NamedTemporaryFile(suffix = '_t.bam', delete=False)
            c = tempfile.NamedTemporaryFile(suffix = '_c.bam', delete=False)
            tname = t.name
            cname = c.name
            t.close()
            c.close()
            bam_sort(args.target_sam, tname, name_sort=True)
            bam_sort(args.contam_sam, cname, name_sort=True)

            # filter to the intermediate output
            f = tempfile.NamedTemporaryFile(suffix = '_f.bam', delete=False)
            fname = f.name
            f.close()
            with pysam.AlignmentFile(tname) as tfile, pysam.AlignmentFile(cname) as cfile:
                filter_file = pysam.AlignmentFile(fname,'wb', template=tfile)
                #contam_file = pysam.AlignmentFile(contam_bam,'wb', template=cfile)
                #unmap_file = pysam.AlignmentFile(unmap_bam,'wb', template=tfile)
                i = 0
                for tread in tfile:
                    cread = cfile.next()
                    assert tread.query_name == cread.query_name
                    if tread.mapping_quality >= min_qual and tread.mapping_quality >= cread.mapping_quality:
                        i += 1
                        filter_file.write(tread)
                    #if tread.mapping_quality < min_qual: # unmapped
                    #    unmap_file.write(tread) # output mapping for target genome
                    #elif tread.mapping_quality < cread.mapping_quality: # contamination
                    #    contam_file.write(cread) # output mapping for contamination genome
                    #else: #target
                    #    filter_file.write(tread)
                filter_file.close()
                #contam_file.close()
                #unmap_file.close()

            # coordinate-sort output
            bam_sort(fname, filter_bam)
            sys.stderr.write('Writing %d alignments to filtered output file %s.\n'%(i, filter_bam))
            #sys.stdout.write('3 output files: filtered %s, contamination %s, unmapped %s\n'%(filter_bam, contam_bam, unmap_bam))

            #clean-up
            os.unlink(tname)
            os.unlink(cname)
            os.unlink(fname)
    else:
        sys.stderr.write('%s filtered alignment exists. OK!\n' % filter_bam)

if __name__ == '__main__':
    main(parse_command_line_arguments())
    
