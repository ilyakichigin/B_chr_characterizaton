#!/usr/bin/env python


import sys
import argparse
from numpy import std


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Outputs region out of position bed file
					usage - get_region.py pos_bed_file.bed genome_file.size > reg.bed
					"""
					)
	parser.add_argument("pos_bed_file", help="bed file with all positions (.bed)")
	parser.add_argument("sizes_file", 
			    help="file which contains size of all chromosomes (.size)")
	parser.add_argument("-d", "--distance", type=int, default=58000,
		            help="restriction on pairwise distance which will show if position in region or not")
	parser.add_argument("-tr", "--total_reads", type=int, default=70,
		            help="minimum number of reads in whole region needed to consider it a true region")

	return parser.parse_args()


def file_into_list(in_file_name):
	
	# Function which transfers all data from file into list

	master_list = [] # This list will contain data from file
	with open(in_file_name, 'rU') as in_file:
		sys.stderr.write("Processing %s\n" % (in_file_name))
		for line in in_file:
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list)

		return master_list


def separ(master_list, n):

	# Using this function we will separate list with positions into many lists containing elements per chromosome

	n = str(n)
	slave_list = []

	for element in master_list:

		if element[0] == 'chr'+n:
			slave_list.append(element)

	return slave_list


def chromosome(n):

	# This function creates simple list with numbers of all chromosomes, depending on n and adds X to the end

	chr_list = range(1,(n+1))
	chr_list.append('X')
	return chr_list
	
		
def delta(p_list):

	# Here position pairwise distance is calculated

	list_with_dist = []
	l = 0

	for element in p_list:

		r = int(element[1])
		if l > 0:
			list_with_dist.append(r-l)
		l = int(element[2])
	return list_with_dist

def average(some_list, n):

	# Calculates mean value of column n

	mean = 0
	j = 0
	total = len(some_list)
	while j < total:
		if n != -1: # If input list is list of lists
			mean += int(some_list[j][n])
		else:
			mean += (some_list[j])
		j += 1
	mean = float(mean) / total
	#sys.stderr.write("\n%f" % (mean))
	return mean
		

if __name__ == '__main__':
	args = parse_command_line_arguments()

	genome_list = file_into_list(args.sizes_file) # Transfer all data from .sizes file into list
	n_chr = len(args.sizes_file)-1 # Number of chromosome not counting X

	pos_list = file_into_list(args.pos_bed_file)

	list_with_chr = chromosome(n_chr)

	#print "chr#\tstart\tend\t#pos\t#reads\tmean PD\tsd\tmean cov\tsd" 

	for n in list_with_chr:

		# In this cycle all actions will be completed per chromosome

		pos_chr_list = separ(pos_list, n) # Get list with positions in chromosome n	
		delta_list = delta(pos_chr_list) # Get list with pairwise distance in chromosome n

		i = 1 # Counter which shows number of position used
		b1 = 0 # Border 1 or start of region
		b2 = 0 # Broder 2 or end of region
		sum_n = 0 # Sum number of reads in a region

		while i < (len(delta_list)-3):

			"""In this cycle we will move from position to position using i
			(actually we are moving from distance to distance,
			 taking into account 2 distances on the left and 2 distances on the right)
			"""

			i += 1 # Because we take distance between 2 on the left and two on the right we should start at +2 position

			k = delta_list[i] # Starting distance
			l = delta_list[i-1] # Distance on the left
			l2 = delta_list[i-2] # Distance on the far left
			r = delta_list[i+1] # Distance on the right
			r2 = delta_list[i+2] # Distance on the far right

			start = pos_chr_list[i-1][1] # Start of region relatively to position we are in (left side relatively to distance k)
			end = pos_chr_list[i+2][2] # End of region relatively to position we are in

	 
			if k < args.distance and \
			   r < args.distance and \
			   r2 < args.distance and \
			   l < args.distance and \
			   l2 > args.distance and \
			   (k+r+r2+l) < (args.distance-0.1*args.distance) and \
			   b1 == 0 and \
			   b2 == 0:

				"""Finding the start of region - all distances should be small except far left one (l2),
				here and below (k+r+(r2)+l+(l2)) parameter will be used as it shows that sum distances in region arent big
				"""

				b1 = start
				sum_pos = 1
				sum_n = 0 # New region => need to nullify sum number of reads in a region
				reads_list = []
				pairwise_dist_list = []
				sum_n += int(pos_chr_list[i-1][3]) + int(pos_chr_list[i-2][3]) # Add number of reads from position we are in and lefter
				reads_list.append(int(pos_chr_list[i-2][3]))
				reads_list.append(int(pos_chr_list[i-1][3]))
				pairwise_dist_list.append(l)
				pairwise_dist_list.append(k)
				

			elif k < args.distance and \
			     r < args.distance and \
			     r2 < args.distance and \
			     l < args.distance and \
			     l2 < args.distance and \
			     (k+r+r2+l+l2) < (args.distance+0.2*args.distance) and \
			     b1 != 0 and \
			     b2 == 0:

				# Elongation of region

				sum_n += int(pos_chr_list[i-1][3])
				reads_list.append(int(pos_chr_list[i-1][3]))
				pairwise_dist_list.append(k)
				sum_pos += 1

			elif k < args.distance and \
			     r < args.distance and \
			     r2 > args.distance and \
			     l < args.distance and \
			     l2 < args.distance and \
			     (k+r+l+l2) < (args.distance-0.1*args.distance) and \
			     b1 != 0 and \
			     b2 == 0:
				sum_n += int(pos_chr_list[i-1][3]) + int(pos_chr_list[i][3]) + int(pos_chr_list[i+1][3])
				reads_list.append(int(pos_chr_list[i-1][3]))
				reads_list.append(int(pos_chr_list[i][3]))
				reads_list.append(int(pos_chr_list[i+1][3]))
				pairwise_dist_list.append(k)
				pairwise_dist_list.append(r)
				if sum_n > args.total_reads:

					"""End of region -  all distances should be small except far right one (r2), region must have start (b1 != 0),
					sum number of reads in a region must not be small and sum distances in region can not be big
					"""

					b2 = end
					sum_pos += 3
					mean_pairwise_dist = sum(pairwise_dist_list)/float(len(pairwise_dist_list))
					mean_cov = sum(reads_list)/float(len(reads_list))
					#print "chr%s\t%s \t%s \t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f" % (n, b1, b2, sum_pos, sum_n, mean_pairwise_dist, std(pairwise_dist_list), mean_cov, std(reads_list))
					print "chr%s\t%s\t%s" % (n, b1, b2)	 
				b1 = 0 # Region ended so we need to nulify all variables
				b2 = 0

	sys.stderr.write("Complete!\n")

