#!/usr/bin/env python

import pysam
import sys
import os.path
import argparse
import subprocess

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes contamination reads by comparing MAPQs to target and cotamination genomes.
                    Input: sam/bam alignments to target and contamination genomes.
                    Output: bam alignments: *filter.bam - specific to target, *contam.bam - specific to contamination, *unmap.bam - unmapped to both.   
                    Pysam python package and samtools are required (tested on v.0.8.1 and v.0.1.19-44428cd, respectively)
                    1) Sort input alignments by read name without cleanup.
                    2) Perform contamination analysis.
                    3) Sort and index resulting alignments with cleanup.
                    4) Remove input files - all read positions are preserved within three generated files
                    """
                    )
    parser.add_argument("target_file",
                        help="sam/bam alignment to target genome")
    parser.add_argument("contam_file",
                        help="sam/bam alignment to contamination genome")
    parser.add_argument("-m", '--min_quality', default=20,
                        help="Minimum quality for filtered file. Default - 20.")
    return parser.parse_args()

   
def sort_by_read_name(filename):
    srt_name = filename[:-4] + '.ns.bam' # ns for name sorted   
    assert os.path.exists(srt_name) == False # input bam not name sorted yet
    sys.stderr.write('samtools view -bS %s | samtools sort -n - > %s\n'%(filename, srt_name))
    p1 = subprocess.Popen(('samtools', 'view', '-bS', filename), stdout=subprocess.PIPE)
    with p1.stdout, open(srt_name, 'w') as outfile:
        p2 = subprocess.Popen(('samtools', 'sort', '-n', '-', srt_name[:-4]), stdin=p1.stdout, stdout=subprocess.PIPE)
    status=[p1.wait(),p2.wait()]
    if p1.returncode != 0 or p1.returncode != 0:
        sys.exit(1)
    os.remove(filename) # remove unsorted file
    sys.stderr.write('rm %s\n'%(filename))
    return srt_name

def sort_index(filename):
    srt_name = filename[:-11] + '.bam' # remove 'unsort.bam'  
    if not os.path.exists(srt_name): # output bam not sorted yet
        sys.stderr.write('samtools sort %s %s; samtools index %s \n'%(filename,srt_name[:-4],srt_name))
        p1 = subprocess.Popen(('samtools', 'sort', filename, srt_name[:-4]), stdout=subprocess.PIPE)
        status=p1.wait()
        if p1.returncode != 0:
            sys.exit(1)
        p2 = subprocess.Popen(('samtools', 'index', srt_name), stdout=subprocess.PIPE)
        status=p2.wait()
        if p2.returncode != 0:
            sys.exit(1)
        os.remove(filename) # remove unsorted file
        sys.stderr.write('rm %s\n'%(filename))
        return srt_name
    else:
        print 'Sorted output alignment %s exists. OK!'%(srt_name)

def compare_mapq(tname, cname, min_qual = 20):
    # inputs: sorted_bam_target.ns.bam, sorted_bam_contam.ns.bam
    # read files - autodetect format, "rb" not specified
    filter_filename = tname[:-7]+'.filter.unsort.bam'   
    contam_filename = cname[:-7]+'.contam.unsort.bam'   
    unmap_filename = tname[:-7]+'.unmap.unsort.bam'   
   
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
    os.remove(tname) # remove input files
    os.remove(cname)
    sys.stderr.write('rm %s %s\n'%(tname, cname))

    return (filter_filename, contam_filename, unmap_filename)   

def main(args):
    args.min_quality = int(args.min_quality)
    args.target_file = sort_by_read_name(args.target_file)
    args.contam_file = sort_by_read_name(args.contam_file)
    outnames = compare_mapq(args.target_file, args.contam_file, args.min_quality)
    for filename in outnames:
        sort_index(filename)

if __name__ == '__main__':
    main(parse_command_line_arguments())
    
