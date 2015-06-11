#!/usr/bin/env python


import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Adds homologies and their scores for each scaffold
					to table from tsv file as last row from homology score file
					"""
					)
	parser.add_argument("tsv_file", help="tsv file")
	parser.add_argument("homology_score_file", help="txt file (obtained by get_hom script) containing scaffolds and their homologies")

	return parser.parse_args()


def addition(master_file_name, second_file_name):
	
	# Adds second column of second file to master file according to scaffolds names

	master_list = [] # This list will contain data from file
	with open(master_file_name, 'rU') as master_file:
		sys.stderr.write("Processing %s\n" % (master_file_name))
		for line in master_file:
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			master_list.append(element_list)
	with open(second_file_name, 'rU') as second_file:
		sys.stderr.write("Processing %s\n" % (second_file_name))
		for line in second_file:
			line = line.strip('\n').strip('\r')
			element_list = line.split('\t')
			i = 0
			while i < len(master_list):
				if element_list[0] == master_list[i][0]:
					master_list[i].append(element_list[1])
				i+=1

		return master_list
	

if __name__ == '__main__':
	args = parse_command_line_arguments()

	tsv_file_name = args.tsv_file
	homology_file_name = args.homology_score_file

	output_list = addition(tsv_file_name, homology_file_name)

	ref_name = homology_file_name.split('.')[1]
	output_list[0].append(ref_name)
	
	for item in output_list:
		print '\t'.join(item)

	sys.stderr.write("Complete!\n")


