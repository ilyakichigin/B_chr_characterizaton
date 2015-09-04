#!/usr/bin/env python

import pysam
import sys
import os.path
import argparse


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes contamination reads by comparing MAPQs to target and cotamination genomes.
                    Input: sam/bam alignments to target and contamination genomes.
                    Output: bam alignments: *filter.bam - specific to target, *contam.bam - specific to contamination, *unmap.bam - unmapped to both.   
                    Pysam python package is required (tested on v.0.8.1)
                    1) (optional) sort input alignments by read name without cleanup.
                    2) Perform contamination analysis.
                    3) sort and index resulting alignments with cleanup.
                    4) Remove input files - all read positions are preserved within three generated files
                    """
                    )
    parser.add_argument("target_file",
                        help="sam/bam alignment to target genome")
    parser.add_argument("contam_file",
                        help="sam/bam alignment to contamination genome")
    parser.add_argument("-m", '--min_quality', default=20,
                        help="Minimum quality for filtered file. Default - 20.")
    parser.add_argument("-a", "--pre_sort_by_name", action="store_true",
                        help="perform preliminary bam sorting by read name. Do not clean up.")
    #parser.add_argument("-p", "--post_sort_index", action="store_true",
    #                    help="perform post sorting by coordinate and indexing of resulting bams. Do clean up.")

    return parser.parse_args()

   
def sort_by_read_name(filename):
    sys.stderr.write('Sorting input alignment by read name: %s\n'%(filename))
    srt_name = filename[:-4] + '.ns.bam' # ns for name sorted   
    assert os.path.exists(srt_name) == False # input bam not name sorted yet
    pysam.sort('-n', filename, srt_name[:-4])
    sys.stderr.write('Output file: %s\n'%(srt_name))
    return srt_name

def sort_index(filename):
    srt_name = filename[:-11] + '.bam' # remove 'unsort.bam'  
    if not os.path.exists(srt_name): # output bam not sorted yet
        sys.stderr.write('Sorting and indexing output alignment: %s\n'%(filename))        
        pysam.sort(filename, srt_name[:-4])
        pysam.index(srt_name)
        sys.stderr.write('Output file: %s\n'%(srt_name))           
        os.remove(filename) # remove unsorted file
        sys.stderr.write('rm %s\n'%(filename))
        return srt_name
    else:
        print 'Sorted output alignment %s exists. OK!'%(srt_name)

def compare_mapq(tname, cname, min_qual = 20, pre_sort_by_name = True):
    # read files - autodetect format, "rb" not specified
    if pre_sort_by_name:
        tst = -7
    else:
        tst = -4
    filter_filename = tname[:tst]+'.filter.unsort.bam'   
    contam_filename = cname[:tst]+'.contam.unsort.bam'   
    unmap_filename = tname[:tst]+'.unmap.unsort.bam'   
   
    with pysam.AlignmentFile(tname) as tfile, pysam.AlignmentFile(cname) as cfile:
        sys.stderr.write('Comparing MAPQs. Target: %s, Contam: %s\n'%(tname,cname))        
        filter_file = pysam.AlignmentFile(filter_filename,'wb', template=tfile)
        contam_file = pysam.AlignmentFile(contam_filename,'wb', template=cfile)
        unmap_file = pysam.AlignmentFile(unmap_filename,'wb', template=tfile)   
        for tread in tfile:
            cread = cfile.next() # gets same line from contamination file
            assert tread.query_name == cread.query_name
            if tread.mapping_quality < min_qual: # unmapped
                unmap_file.write(tread) # output mapping for target genome
            elif tread.mapping_quality < cread.mapping_quality: # contamination
                contam_file.write(cread) # output mapping for contamination genome
            else: #target
                filter_file.write(tread)

        filter_file.close()
        contam_file.close()
        unmap_file.close()

    sys.stderr.write('Output files: %s, %s, %s\n'%(filter_filename, contam_filename, unmap_filename))
    # remove input files
    os.remove(tname)
    os.remove(cname)
    sys.stderr.write('rm %s %s\n'%(tname, cname))

    return (filter_filename, contam_filename, unmap_filename)   

def main(args):
    args.min_quality = int(args.min_quality)

    if args.pre_sort_by_name:
        args.contam_file = sort_by_read_name(args.contam_file)
        args.target_file = sort_by_read_name(args.target_file)

    outnames = compare_mapq(args.target_file, args.contam_file, args.min_quality, args.pre_sort_by_name)
   
    for filename in outnames:
        sort_index(filename)
    
    os.remove(args.contam_file)
    os.remove(args.target_file)

if __name__ == '__main__':
    main(parse_command_line_arguments())
    
