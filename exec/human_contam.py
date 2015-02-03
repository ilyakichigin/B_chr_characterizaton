#!/usr/bin/env python

import pysam
import sys
import os.path
import argparse

#input
#file1 - bam/sam vs target genome sorted by read name
#file2 - bam/sam vs contamination genome sorted by read name

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(prog="compEDps",description=    
                    """
                    Removes contamination reads by comparing MAPQs to target and cotamination genomes.
                    Input: sam/bam alignments to target and contamination genomes.
                    Output: bam alignments: *filter.bam - specific to target, *contam.bam - specific to contamination, *unmap.bam - unmapped to both.   
                    Pysam python package is required (tested on v.0.8.1)
                    1) (optional) sort input alignments by read name without cleanup.
                    2) Perform contamination analysis.
                    3) (optional) sort and index resulting alignments with cleanup.
                    """
                    )
    parser.add_argument("target_file",
                        help="sam/bam alignment to target genome")
    parser.add_argument("contam_file",
                        help="sam/bam alignment to contamination genome. Default - 20.")
    parser.add_argument("-m", '--min-quality', default=20,
                        help="Minimum quality to output in filtered filtered")
    parser.add_argument("-a", "--pre-sort-by-name", action="store_true",
                        help="perform preliminary bam sorting by read name. Do not clean up.")
    parser.add_argument("-p", "--post-sort-index", action="store_true",
                        help="perform post sorting by coordinate and indexing of resulting bams. Do clean up.")

    return parser.parse_args()

   
def sort_by_read_name(filename):
    srt_name = filename[:-4] + '.ns.bam' # ns for name sorted   
    assert os.path.exists(srt_name) == False
    sys.stdout.write('Sorting input alignment by read name: %s\n'%(filename))
    pysam.sort('-n', filename, srt_name[:-4])
    sys.stdout.write('Output file: %s\n'%(srt_name))
    return srt_name

def sort_index(filename):
    srt_name = filename[:-4] + '.cs.bam' # cs for coordinate sorted   
    assert os.path.exists(srt_name) == False
    sys.stdout.write('Sorting and indexing output alignment: %s\n'%(filename))
    pysam.sort(filename, srt_name[:-4])
    pysam.index(srt_name)
    os.remove(filename)
    sys.stdout.write('Output file: %s\n'%(srt_name))   
    return srt_name


def compare_mapq(tname, cname, min_qual = 20):
    # read files - autodetect format, "rb" not specified
    if args.pre_sort_by_name:
        tst = -7
    else:
        tst = -4
    filter_filename = tname[:tst]+'.filter.bam'   
    contam_filename = cname[:tst]+'.contam.bam'   
    unmap_filename = tname[:tst]+'.unmap.bam'   
   
    with pysam.AlignmentFile(tname) as tfile, pysam.AlignmentFile(cname) as cfile:
        filter_file = pysam.AlignmentFile(filter_filename,'wb', template=tfile)
        contam_file = pysam.AlignmentFile(contam_filename,'wb', template=cfile)
        unmap_file = pysam.AlignmentFile(unmap_filename,'wb', template=tfile)   
        sys.stdout.write('Comparing MAPQs. Target: %s, Contam: %s\n'%(tname,cname))        
        for tread in tfile:
            cread = cfile.next() # gets same line from contamination file
            assert tread.query_name == cread.query_name
            if (tread.mapping_quality < min_qual) and (cread.mapping_quality < min_qual): # unmapped
                unmap_file.write(tread) # output mapping for target genome
            elif tread.mapping_quality < cread.mapping_quality: # contamination
                contam_file.write(cread) # output mapping for contamination genome
            else: #target
                filter_file.write(tread)

        filter_file.close()
        contam_file.close()
        unmap_file.close()

    sys.stdout.write('Output files: %s, %s, %s\n'%(filter_filename, contam_filename, unmap_filename))

    return (filter_filename, contam_filename, unmap_filename)   

if __name__ == '__main__':
    args = parse_command_line_arguments()

    if args.pre_sort_by_name:
        args.contam_file = sort_by_read_name(args.contam_file)
        args.target_file = sort_by_read_name(args.target_file)

    outnames = compare_mapq(args.target_file, args.contam_file, args.min_quality)
   
    if args.post_sort_index:
        for filename in outnames:
           sort_index(filename)
