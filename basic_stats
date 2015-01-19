#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Creates reads, positions, reads/positions in region and out of region
					bed files from input reads bam file, region bed file and genome file
					after that calculates and sends statistic to stdout.
					"""
					)
	parser.add_argument("bam_file", help="bam file with reads (.bam)")
	parser.add_argument("region_file", help="bed file containing regions (.reg.bed)")
	parser.add_argument("genome_file", help="genome file which contains size of all chromosomes and total genome size at the end (.genome)")

	return parser.parse_args()


def creating_files(arg_list, c):

	# Function which is used to create .bed files through bedtools commands

	proc_list = ['bedtools'] + arg_list
	process = subprocess.Popen(proc_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(output, err) = process.communicate()
	log = 'a' # log shows if we need to create log file anew or add something to existent one

	if err == "": # If statement that shows if bedtools did all with no errors and stderr output is empty or not 
		status="completed"
	else:
		status="attempted"

	if c == 0: # Conditions which control name of created file and log message depending on number of command - c
		ending = ".reads.bed"
		err_message = "Bam to bed command %s\n" % (status)
		log = 'w' # First log file should be created anew
	elif c == 1:
		ending = ".pos.bed"
		err_message = "Position creation %s\n" % (status)
	elif c == 2:
		ending = ".readsin.bed"
		err_message = "Finding reads in region %s\n" % (status)
	elif c == 3:
		ending = ".posin.bed"
		err_message = "Position creation in region %s\n" % (status)
	elif c == 4:
		ending = ".readsout.bed"
		err_message = "Finding reads out of region %s\n" % (status)
	elif c == 5:
		ending = ".posout.bed"
		err_message = "Position creation out of region %s\n" % (status)

	out_file_name = in_file_name[:-4] + ending # This section represent creation of files
	out_file = open(out_file_name, 'w')
	out_file.write(output)
	out_file.close()
	file_list.append(out_file_name) # List with names of all created files

	err_file_name = in_file_name[:-4] + ".log" # This section represent creation of log file
	err_file = open(err_file_name, log)
	err_file.write(err + err_message)
	err_file.close()


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

	
def calculating_bp_region(reg):

	# This function calculates bp size of region from region .bed file

	s = 0
	with open(reg, 'rU') as in_file:
		sys.stderr.write("Processing %s\n" % (reg))
		for line in in_file:
			line = line.strip('\n').strip('\r')
			element_list=line.split('\t')
			s += (int(element_list[2]) - int(element_list[1]))
		
		bp_region = float(s)
		return bp_region


def get_sum(some_list):

	# Small function to calculate sum bp in list which was acquired from .bed file

	s = 0
	for element in some_list[1:]: # Note that first element of list ignored because it contains name of file
		s += (int(element[2]) - int(element[1]))

	return s


def get_stats(in_file_name, bp):

	# Function that calculates number, sum bp and density of reads/positions and outputs it into stdout

	with open(in_file_name, 'rU') as in_file:
		sys.stderr.write("Processing %s\n" % (in_file_name))
		master_list = [] # Input .bed file will be loaded to this list
		name = in_file_name.split('.')[-2] # This name will be used to mark master_list, it will be like ending before .bed in creating_files function
		master_list.append(name)

		for line in in_file: # Cycle to convert all lines of file to master_list, making list of lists
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list)

	number = len(master_list) - 1 # Number of reads or positions, 1 is substracted because master_list[0] contains name of file
	s = get_sum(master_list) # Sum bp of reads or positions

	print "Number of %s\t%5d" % (master_list[0], number)
	print "Sum bp %s\t%5d" % (master_list[0], s)
	print "%s density\t%5.8f\n" % (master_list[0], float(number)/bp)

	return master_list


def in_out (in_list, out_list, bp_in, bp_out):

	# Function which compares reads or positions stats in region to out of region

	s_in = get_sum(in_list)
	s_out = get_sum(out_list)
	number_in = float(len(in_list)-1)
	number_out = float(len(out_list)-1)
	
	print "Number of %s/%s of region\t%5.8f" % (in_list[0], out_list[0], number_in/number_out)
	print "Sum bp of %s/%s of region\t%5.8f" % (in_list[0], out_list[0], float(s_in)/s_out)
	print "Density of %s/%s of region\t%5.8f\n" % (in_list[0], out_list[0], (number_in/bp_in)/(number_out/bp_out))

	
def cov_sd (r_list, p_list):

	#This function calculates coverage and it's standard deviation

	cov = (float((len(r_list)-1))) / (float((len(p_list)-1)))

	s = 0
	for element in p_list[1:]:
		s += (float(element[3]) - cov) ** 2		
	sd = (s / (len(p_list)-1)) ** 0.5

	print "Coverage %s/%s\t%5.8f" % (r_list[0], p_list[0], cov)
	print "Standard deviation\t%5.8f\n" % (sd)
	

if __name__ == '__main__':
	args = parse_command_line_arguments()

	in_file_name = args.bam_file
	region_file = args.region_file
	genome_file = args.genome_file

	assert in_file_name.endswith('.bam')
	assert region_file.endswith(".reg.bed")
	assert genome_file.endswith(".genome")

	file_list = [] # List which will contain all created files names
	i = 0 # This variable will show number of command (0 - first, 1 - second and so on)

	command_list =	[
			['bamtobed', '-i', in_file_name],
			['merge', '-n', '-i', in_file_name[:-4]+".reads.bed"],
			['intersect', '-a', in_file_name[:-4]+".reads.bed", '-b', region_file],
			['merge', '-n', '-i', in_file_name[:-4]+".readsin.bed"],
			['intersect', '-v', '-a', in_file_name[:-4]+".reads.bed", '-b', region_file],
			['merge', '-n', '-i', in_file_name[:-4]+".readsout.bed"],
			]

	while i < 6: # Through this cycle 6 .bed files will be created - reads, positions, reads in region, positions in region, reads out of region and positions out of region
		creating_files(command_list[i], i)
		i += 1

	bp_genome = float(file_into_list(genome_file)[-1][1]) # Size of genome is taken from genome file as last line
	bp_region = calculating_bp_region(region_file) # Use function to calculate size of region

	reads_list = get_stats(file_list[0], bp_genome) # With this section basic statistic will be outputed and data in bed files will be copied into lists to use it later
	pos_list = get_stats(file_list[1], bp_genome)
	readsin_list = get_stats(file_list[2], bp_region)
	posin_list = get_stats(file_list[3], bp_region)
	readsout_list = get_stats(file_list[4], (bp_genome-bp_region))
	posout_list = get_stats(file_list[5], (bp_genome-bp_region))

	in_out(readsin_list, readsout_list, bp_region, (bp_genome-bp_region)) # Functions to compare reads/positions stats in region and out of region using previously created lists
	in_out(posin_list, posout_list, bp_region, (bp_genome-bp_region))

	cov_sd(reads_list,pos_list) # Functions to calculate coverage and it's standard deviation using previously created lists
	cov_sd(readsin_list,posin_list)
	cov_sd(readsout_list,posout_list)

	sys.stderr.write("Completed\n")
