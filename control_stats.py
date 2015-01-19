#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Creates reads and positions bed files from input reads bam file
					and genome file after that creates bed files for each
					individual chromosome and calculates and outputs their statistics
					in form of a table .
					"""
					)
	parser.add_argument("bam_file", help="bam file with reads (.bam)")
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


def separ(master_list, n):

	# That function separates list with positions/reads into many lists containing elements per chromosome

	n = str(n)
	slave_list = []
	ending="." + master_list[0] + '.chr' + n + '.bed' # Ending that will be added to end of file, master_list[0] holds either reads or pos

	for element in master_list[1:]:

		if element[0] == 'chr'+n:
			slave_list.append(element)

	out_file_name = in_file_name[:-4] + ending
	out_file = open(out_file_name,'w')
	sys.stderr.write("Creating %s\n" % (out_file_name))
	output = ''

	for element in slave_list:	# Save chromosome data in file
		output += '\t'.join(element) + "\n"

	out_file.write(output)
	out_file.close()

	return slave_list


def chromosome(n):

	# This function creates simple list with numbers of all chromosomes, depending on n and adds X to the end

	chr_list = range(1,n)
	chr_list.append('X')
	return chr_list

				
def get_sum(some_list):

	# Small function to calculate sum bp in list which was acquired from .bed file

	s = 0
	for element in some_list:
		s += (int(element[2]) - int(element[1]))

	return s


def get_list(in_file_name):

	# Function that creates list out of .bed file

	with open(in_file_name, 'rU') as in_file:
		sys.stderr.write("Processing %s\n" % (in_file_name))
		master_list = [] # Input .bed file will be loaded to this list
		name = in_file_name.split('.')[-2] # This name will be used to mark master_list, it will be like ending before .bed in creating_files function
		master_list.append(name)

		for line in in_file: # Cycle to convert all lines of file to master_list, making list of lists
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list)

	return master_list


def get_stats(r_list, p_list, n):

	# Function that calculates statistics per chromosome

	r_number = len(r_list)	# Number of reads/positions
	p_number = len(p_list)
	r_s = get_sum(r_list) # Sum bp of reads/positions
	p_s = get_sum(p_list)
	final_list.append('chr' + str(n) + '\t') # Add new row to table-list which will hold chromosome statistic and will be outputed in the end as table, n shows number of row/chromosome
	if n == 'X':
		n = len(genome_list) - 2	# If n = X it is not the number, so need to assign to n X's position which is last before total

	final_list[n] += '%d\t' % (r_number) # Add chromosome statistics to table-list
	final_list[n] += '%d\t' % (r_s)
	final_list[n] += '%.6f\t' % (r_number/float(genome_list[int(n)-1][1]))
	final_list[n] += '%d\t' % (p_number)
	final_list[n] += '%d\t' % (p_s)
	final_list[n] += '%.6f\t' % (p_number/float(genome_list[int(n)-1][1]))

	return (r_list, p_list)


	
def delta_sd(r_list, p_list, n):

	# Here mean difference between reads/positions and it's standard deviation is calculated

	delta_list = []
	mean_delta = 0
	s = 0
	l = 0
	if n == 'X':
		n = len(genome_list) - 2

	for element in r_list:

		r = int(element[1])
		if l > 0:
			delta_list.append(r-l)
		l = int(element[2])

	for delta in delta_list:
		mean_delta += delta
	mean_delta = float(mean_delta)/len(delta_list)

	for delta in delta_list:
		s += (float(delta) - mean_delta) ** 2
	sd = (s / len(delta_list)) ** 0.5

	final_list[n] += '%.1f\t' % (mean_delta) # Add chromosome statistics to table-list
	final_list[n] += '%.1f\t' % (sd)

	delta_list = []
	mean_delta = 0
	s = 0
	l = 0

	for element in p_list:

		r = int(element[1])
		if l > 0:
			delta_list.append(r-l)
		l = int(element[2])

	for delta in delta_list:
		mean_delta += delta
	mean_delta = float(mean_delta)/len(delta_list)

	for delta in delta_list:
		s += (float(delta) - mean_delta) ** 2
	sd = (s / len(delta_list)) ** 0.5

	final_list[n] += '%.1f\t' % (mean_delta) # Add chromosome statistics to table-list
	final_list[n] += '%.1f\t' % (sd)


def cov_sd (r_list, p_list, n):

	#This function calculates coverage and it's standard deviation

	if n == 'X':
		n = len(genome_list) - 2
	s = 0
	for element in p_list[1:]:
		s += float(element[3])
	cov = s/(len(p_list)-1)

	s = 0
	for element in p_list[1:]:
		s += (float(element[3]) - cov) ** 2		
	sd = (s / (len(p_list)-1)) ** 0.5

	final_list[n] += '%.3f\t' % (cov) # Add chromosome statistics to table-list
	final_list[n] += '%.3f\t' % (sd)


if __name__ == '__main__':
	args = parse_command_line_arguments()

	in_file_name = args.bam_file
	genome_file = args.genome_file

	assert in_file_name.endswith('.bam')
	assert genome_file.endswith('.genome')

	file_list = [] # List which will contain all created files names
	i = 0 # This variable will show number of command (0 - first, 1 - second)

	''' Table-list which will hold statistics of each individual chromosome on each row,
	first row containing names of values is created at start'''
	final_list = ['\tN reads\tRead bp\tRead dens\tN pos\tPos bp\tPos dens\tRead mean delta\tStandard dev\tPos mean delta\tStandard dev\tCoverage\tStandard dev\t']

	genome_list = file_into_list(genome_file) # Transfer all data from .genome file into list
	n_chr = len(genome_list)-2 # Number of chromosome not counting X and total

	command_list = 	[
			['bamtobed', '-i', in_file_name],
			['merge', '-n', '-i', in_file_name[:-4]+".reads.bed"]
			]

	while i < 2: # Through this cycle 2 .bed files will be created - reads and positions
		creating_files(command_list[i], i)
		i += 1

	reads_list = get_list(file_list[0])
	pos_list = get_list(file_list[1])

	list_with_chr = chromosome(n_chr)
	
	for n in list_with_chr: # In this cycle lists will be parsed through series of functions to get all needed stats

		(_in_r_list, _in_p_list) = get_stats(separ(reads_list, n), separ(pos_list, n), n)
		delta_sd(_in_r_list, _in_p_list, n)

		cov_sd(get_list(in_file_name[:-4] + '.reads.chr' + str(n) + '.bed'), get_list(in_file_name[:-4] + '.pos.chr' + str(n) + '.bed'), n)

	table = '\n'.join(final_list)
	print table


	sys.stderr.write("Completed\n")
