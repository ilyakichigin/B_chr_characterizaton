#!/usr/bin/env python


import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Obtains sizes of random ACA6 scaffolds
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
			if element_list[0] == 'chr6':
				mediate_list =[]
				mediate_list.append(element_list[3])
				mediate_list.append(str(int(element_list[2])-int(element_list[1])))
				master_list.append(mediate_list)
		return master_list

if __name__ == '__main__':
	args = parse_command_line_arguments()

	output_list = file_into_list(args.input_file)

	for item in output_list:
		print '\t'.join(item)

	sys.stderr.write("Complete!\n")
