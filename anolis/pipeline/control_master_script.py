#!/usr/bin/env python

import sys
import argparse
import os

def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=    
		"""
		Master script for generating random scaffolds on chromosome 6 of A. carolinensis
		1)Generating bed files with random scaffolds out of scaffold size file
		2)Adding file generated to bed with whole genome excluding chromosome 6
		3)Sorting new bedfile		
		4)Bedtools getfasta to obtain .fa file with new scaffolds		
		5)Bowtie2 indexing of new fasta file		
		6)Bowtie2 align to target and contaminate genomes
		7)Filtering sams obtained
		8)Bedtools commands to get reads file
		9)Getting file with new scaffolds sizes
		10)Cleaning up
		"""
		)
	parser.add_argument("scaffolds_sizes_file", help="file with scaffolds sizes")

	parser.add_argument("-i", "--number_iterations", default="10", help="number of iterations of this script")

	parser.add_argument("--path_to_bedtools", default="bedtools", help="path to bedtools binary")

	parser.add_argument("--path_to_bowtie2", default="bowtie2", help="path to bowtie2 binary")

	parser.add_argument("-p", "--proc_bowtie2", default="1", help="number of processors allocated for bowtie2. Default - 1.")

	parser.add_argument("--path_to_anoCar2.fa", default="anoCar2", help="path to reference anoCar2.fa")

	return parser.parse_args()


if __name__ == '__main__':
	args = parse_command_line_arguments()

	for i in range(1, int(args.number_iterations)+1):

		command_list = [
			('exec/random_scaffolds.py ' + args.scaffolds_sizes_file + ' > random_scaffolds.bed'),
			('cat random_scaffolds.bed anoCar2_no_chr6.bed > anoCar2_new_chr6.bed'),
			(args.path_to_bedtools + ' sort -i anoCar2_new_chr6.bed > anoCar2_new_chr6.srt.bed'),
			(args.path_to_bedtools + ' getfasta -fi ' + args.path_to_anoCar2.fa + ' -bed anoCar2_new_chr6.srt.bed -fo anoCar2_new_chr6.fa -name'),
			(args.path_to_bowtie2 + '-build anoCar2_new_chr6.fa anoCar2_new_chr6'),
			(args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x anoCar2_new_chr6 -1 ACAF.F.ca.fastq -2 ACAF.R.ca.fastq -S ACAF_control.anoCar2.sam'),
			(args.path_to_bowtie2 + ' -p ' + args.proc_bowtie2 + ' --local -x genome -1 ACAF.F.ca.fastq -2 ACAF.R.ca.fastq -S ACAF_control.hg19.sam'),
			('exec/contam_filter.py ACAF_control.anoCar2.sam ACAF_control.hg19.sam -a'),
			(args.path_to_bedtools + ' bamtobed -i ACAF_control.anoCar2.filter.bam > ACAF_control.' + str(i) + '.reads.bed'),
			('exec/get_scaffold_size.py anoCar2_new_chr6.srt.bed > random_scaffolds.sizes'),
			('cat random_scaffolds.sizes anoCar2_no_chr6.sizes > anoCar2_rs.sizes'),
			('exec/sort_scaffolds.py anoCar2_rs.sizes > anoCar2.control.' + str(i) +'.sizes'),
			('rm random_scaffolds.bed anoCar2_new_chr6.bed anoCar2_new_chr6.srt.bed anoCar2_new_chr6.fa ACAF_control.anoCar2.sam anoCar2*.bt2  ACAF_control.hg19.sam *.bam random_scaffolds.sizes anoCar2_rs.sizes')
			]
	    
		# run
		for command in command_list:
			sys.stderr.write(command+'\n') 
			os.system(command)
    
	sys.stderr.write("Complete!\n")
