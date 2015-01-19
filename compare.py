#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Gets difference between 2 positions files and
					compares number and sum length of positions
					"""
					)
	parser.add_argument("pos1_file", help="main bed file with positions (.bed)")
	parser.add_argument("pos2_file", help="secondary bed file with positions (.bed)")

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
		name = in_file1_name[:-4] + ".intersectwith." + in_file2_name.split('.')[0] + "." + in_file2_name.split('.')[-2] + ".bed"
		err_message = "Intersect %s\n" % (status)
		log = 'w' # First log file should be created anew
	elif c == 1:
		name = in_file1_name[:-4] + ".oppositewith." + in_file2_name.split('.')[0] + "." + in_file2_name.split('.')[-2] + ".bed"
		err_message = "Opposite intersect #1 %s\n" % (status)
	elif c == 2:
		name = in_file2_name[:-4] +  ".oppositewith." + in_file1_name.split('.')[0] + "." + in_file1_name.split('.')[-2] + ".bed"
		err_message = "Opposite intersect #2 %s\n" % (status)

	out_file_name = name # This section represent creation of files
	out_file = open(out_file_name, 'w')
	out_file.write(output)
	out_file.close()
	file_list.append(out_file_name) # List with names of all created files

	err_file_name = in_file1_name[:-4] + in_file2_name.split('.')[0] + ".log" # This section represent creation of log file
	err_file = open(err_file_name, log)
	err_file.write(err + err_message)
	err_file.close()


def get_sum(some_list):

	# Small function to calculate sum bp in list which was acquired from .bed file

	s = 0
	for element in some_list: 
		s += (int(element[2]) - int(element[1]))

	return s


def get_stats(in_file_name):

	# Function that calculates number, sum bp of positions

	with open(in_file_name, 'rU') as in_file:
		#sys.stderr.write("Processing %s\n" % (in_file_name))
		master_list = [] # Input .bed file will be loaded to this list

		for line in in_file: # Cycle to convert all lines of file to master_list, making list of lists
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list)

	number = len(master_list) # Number of positions
	s = get_sum(master_list) # Sum bp of positions
	
	return (number, s)


if __name__ == '__main__':
	args = parse_command_line_arguments()

	in_file1_name = args.pos1_file
	in_file2_name = args.pos2_file

	assert in_file1_name.endswith('.bed')
	assert in_file2_name.endswith(".bed")



	file_list = [] # List which will contain all created files names
	i = 0 # This variable will show number of command (0 - first, 1 - second and so on)

	command_list =	[
			['intersect', '-a', in_file1_name, '-b', in_file2_name],
			['intersect', '-v', '-a', in_file1_name, '-b', in_file2_name],
			['intersect', '-v', '-a', in_file2_name, '-b', in_file1_name]
			]

	while i < 3: # Through this cycle 3 .bed files will be created
		creating_files(command_list[i], i)
		i += 1

	(number_same, s_same) = get_stats(file_list[0])
	(number_different1, s_different1) = get_stats(file_list[1])
	(number_different2, s_different2) = get_stats(file_list[2])

	print "Number of positions which are both in %s and in %s\t%5d" % (in_file1_name.split('.')[0], in_file2_name.split('.')[0], number_same)
	print "Sum length of positions which are both in %s and in %s\t%5d\n" % (in_file1_name.split('.')[0], in_file2_name.split('.')[0], s_same)

	print "Number of positions which are in %s but not in %s\t%5d" % (in_file1_name.split('.')[0], in_file2_name.split('.')[0], number_different1)
	print "Sum length of positions which are in %s but not in %s\t%5d\n" % (in_file1_name.split('.')[0], in_file2_name.split('.')[0], s_different1)

	print "Number of positions which are in %s but not in %s\t%5d" % (in_file2_name.split('.')[0], in_file1_name.split('.')[0], number_different2)
	print "Sum length of positions which are in %s but not in %s\t%5d\n" % (in_file2_name.split('.')[0], in_file1_name.split('.')[0], s_different2)

	print "Similarity between number of positions relatively to %s\t%5.2f%%" % (in_file1_name.split('.')[0], (1 - (float(number_different1)/float(number_same+number_different1))) * 100)
	print "Similarity between sum length of positions relatively to %s\t%5.2f%%\n" % (in_file1_name.split('.')[0], ((1 - float(s_different1)/float(s_same+s_different1))) * 100)

	print "Similarity between number of positions relatively to %s\t%5.2f%%" % (in_file2_name.split('.')[0], (1 - (float(number_different2)/float(number_same+number_different2))) * 100)
	print "Similarity between sum length of positions relatively to %s\t%5.2f%%\n" % (in_file2_name.split('.')[0], ((1 - float(s_different2)/float(s_same+s_different2))) * 100)
