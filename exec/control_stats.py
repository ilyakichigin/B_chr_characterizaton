#!/usr/bin/env python


import subprocess
import sys
import argparse
import numpy


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    With bam file and chromosome sizes file (.sizes) as inputs create reads and positions bed file and print statistics for each chromosome. 
                    Requires: bedtools (tested on v.2.17.0), numpy (tested on v.1.8.2) 
                    """
                    )
    parser.add_argument("bam_file", help="input bam file")
    parser.add_argument("sizes_file", help="tab-separated file with sizes of chromosomes (for example, from UCSC genome database chromInfo.txt) #and total genome size at the end# (.sizes). Note that only chromosomes listed in this file will be processed.")
    parser.add_argument("--path_to_bedtools", default="bedtools",help="path to bedtools")

    return parser.parse_args()
    
def parse_sizes(filename):
    
    # Parse chromosome size file. Returns list: [(chromosome, size)]. 
    # remove genome size

    master_data = [] 
    with open(filename, 'rU') as in_file:
        sys.stderr.write("Processing chromosome sizes file %s\n" % (filename))

        for line in in_file:
            element_list = line.split('\t')
            master_data.append((element_list[0], element_list[1]))
            
    return master_data

def run_bedtools(bam_file, path_to_bedtools):

    # Using bedtools generate .bed files with reads and positions. Returns output file names.
    
    file_list = []

    read_bed_name = bam_file[:-4]+'.reads.bed'
    pos_bed_name = bam_file[:-4]+'.pos.bed'   
    
    command_list = [ # [command string, output file name]
            [path_to_bedtools + ' bamtobed -i ' + bam_file, read_bed_name],
            [path_to_bedtools + ' merge -n -i ' + read_bed_name, pos_bed_name]
            ]    

    for command in command_list:
        sys.stderr.write(command[0]+'\n')
        process = subprocess.Popen(command[0].split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output, err) = process.communicate()
        
        sys.stderr.write(err)
    
        with open(command[1],'w') as out_file:
            sys.stderr.write('Writing output to file '+command[1]+'\n')
            out_file.write(output)
        file_list.append(command[1]) # List with names of created files
        
    return file_list


def get_list(in_file_name):

    # Reads file into list of lists. 

    with open(in_file_name, 'rU') as in_file:
        sys.stderr.write("Reading %s to list.\n" % (in_file_name))
        master_list = [] 
        for line in in_file:
            element_list = line.split()
            master_list.append(element_list)

    return master_list


def get_basic_stats(bed_list, chrom_size):

    # Calculate basic statistics for list from .bed file
    
    # Number of reads/positions
    bed_num = len(bed_list) 
    # Sum bp of reads/positions
    bed_sum = sum([(int(x[2]) - int(x[1])) for x in bed_list]) 
    # Density
    bed_dens = bed_num/float(chrom_size)
    
    return (bed_num, bed_sum, bed_dens) 

    
def calc_delta(bed_list):

    # Calculate mean and sd distance between reads/positions for list from .bed file

    # generate list of distances between neighbors    
    delta_list = []
    l = 0
    for element in bed_list:
        r = int(element[1])
        if l > 0:
            delta_list.append(r-l)
        l = int(element[2])

    # calculate mean and standard deviation
    mean_delta = numpy.mean(delta_list)
    sd_delta = numpy.std(delta_list)

    return (mean_delta, sd_delta)


def calc_cov(bed_list):

    # Calculate mean and sd coverage of positions for list from .bed file

    cov_list = [int(element[3]) for element in bed_list]
    mean_cov = numpy.mean(cov_list)
    sd_cov = numpy.std(cov_list)   
    
    return (mean_cov, sd_cov)

if __name__ == '__main__':
    args = parse_command_line_arguments()

    assert args.bam_file.endswith('.bam')
    assert args.sizes_file.endswith('.sizes')
    
    # Generate bed files
    wg_bed_files = run_bedtools(args.bam_file, args.path_to_bedtools)                
    read_list = get_list(wg_bed_files[0])
    pos_list = get_list(wg_bed_files[1])
    # Parse size file
    size_data = parse_sizes(args.sizes_file)    
    
    # Calculate and print stats
    header = ['chrom','reads','reads.bp','reads.dens','pos','pos.bp','pos.dens','read.pd.mean','read.pd.sd',
'pos.pd.mean','pos.pd.sd','pos.cov.mean','pos.cov.sd']
    print '\t'.join(header)
    for (chrom, chrom_size) in size_data:
        # subset reads and positions from chromosome
        chrom_r_list = [read for read in read_list if read[0] == chrom]
        chrom_p_list = [pos for pos in pos_list if pos[0] == chrom]
        # calculate stats
        chrom_stats = [chrom]
        chrom_stats.extend( get_basic_stats(chrom_r_list, chrom_size) ) # read basic stats
        chrom_stats.extend( get_basic_stats(chrom_p_list, chrom_size) ) # position basic stats
        chrom_stats.extend( calc_delta(chrom_r_list) ) # read pairwise distance (pd) stats - pretty useless
        chrom_stats.extend( calc_delta(chrom_p_list) ) # position distance stats
        chrom_stats.extend( calc_cov(chrom_p_list) ) # position coverage stats
        chrom_stats = [str(x) for x in chrom_stats]
        print '\t'.join(chrom_stats)
        
    sys.stderr.write("Complete!\n")
