#!/usr/bin/env python


import sys
import argparse
import random


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Creates random scaffolds for ACA chr 6
					"""
					)
	parser.add_argument("input_file", help="scaffold sizes file")

	return parser.parse_args()


def file_into_list(file_name):

	master_list = [] # This list will contain data from file
	with open(file_name, 'rU') as master_file:
		sys.stderr.write("Processing %s\n" % (file_name))
		for line in master_file:
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list[1])
		return master_list


def create_scaffolds(input_list):

	i = 0
	n = 0
	master_list = []
	while i < 80741955: # 80741955 is a size of ACA6
		j = int(random.choice(input_list))
		start = str(i + 1)
  		end = str(j + i)
		if i+j > 80741955:
			end = str(80741955)
		element_list = []
		element_list.append('chr6')
  	element_list.append(start)
		element_list.append(end)
		element_list.append('control_%d' % (n))
		master_list.append(element_list)
		i += j
		n += 1
	return master_list	

if __name__ == '__main__':
	args = parse_command_line_arguments()

	mediate_list = file_into_list(args.input_file)
	output_list = create_scaffolds(mediate_list)

	for item in output_list:
		print '\t'.join(item)
