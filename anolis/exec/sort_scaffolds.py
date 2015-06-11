#!/usr/bin/env python


import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Sort scaffold sizes file descendant according to size
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
			if len(element_list) > 1:
				element_list[1] = int(element_list[1])
				master_list.append(element_list)
		return master_list

if __name__ == '__main__':
	args = parse_command_line_arguments()

	mediate_list = file_into_list(args.input_file)
	output_list = sorted(mediate_list, key=lambda scaffold: scaffold[1], reverse=True)

	for item in output_list:
		item[1] = str(item[1])
		print '\t'.join(item)

	sys.stderr.write("Complete!\n")
