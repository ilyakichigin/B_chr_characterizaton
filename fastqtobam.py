#!/usr/bin/env python


import subprocess
import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Complete series of commands to create and sort bam files
					out of fastq files, reads in direct and reverse files
					will be renamed to show their directivity, paired end
					mapping will be done to the reference genome and human genome (hg19),
					log file with stderr of programms used will be created,
					should be ran in a folder with indexed refernce genomes
					"""
					)
	parser.add_argument("fastq_R1_file", help="fastq file with forward reads (.fastq)")

	parser.add_argument("fastq_R2_file", help="fastq file with reverse reads (.fastq)")

	parser.add_argument("name", help="name of sample (can be anything, in our case CFA12, VVUB and so on)")

	parser.add_argument("name_ref", help="name of reference genome (canFam3, bosTau7)")

	return parser.parse_args()
   

def creating_files(arg_list, c):

	# Function which is used to run programms

	process = subprocess.Popen(arg_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(output, err) = process.communicate()
	log = 'a' # log shows if we need to create log file anew or add something to existent one

	if err == "": # If statement that shows if corresponding programm did all with no messages and stderr output is empty or not 
		status="completed with no messages"
	else:
		status="is done with some messages"

	if c == 0: # c shows the number of process and is needed to contrtol names of output files
		ending = ".R1.ca.fastq"
		err_message = "Cut adapt #1 %s\n" % (status)
		log = 'w' # First log file should be created anew
	elif c == 1:
		ending = ".R2.ca.fastq"
		err_message = "Cut adapt #2 %s\n" % (status)
	elif c == 2:
		ending = ".ca." + name_ref + ".pe.sam"
		err_message = "Alignment %s\n" % (status)
	elif c == 3:
		ending = ".ca.hg19.pe.sam"
		err_message = "Alignment to hg19 %s\n" % (status)		
	elif c == 4:
		ending = ".ca." + name_ref + ".pe.bam"
		err_message = "Sam to bam %s\n" % (status)
	elif c == 5:
		ending = ""
		err_message = "Sorting bam %s\n" % (status)		
	elif c == 6:
		ending = ".ca.hg19.pe.bam"
		err_message = "Sam to bam %s\n" % (status)
	elif c == 7:
		ending = ""
		err_message = "Sorting bam %s\n" % (status)


	if c != 5 and c != 7:
		out_file_name = name + ending # This section represent creation of files (on step 5 and 7 files are created automatically)
		out_file = open(out_file_name, 'w')
		out_file.write(output)
		out_file.close()

	err_file_name = name + ".log" # This section represent creation of log file
	err_file = open(err_file_name, log)
	err_file.write(err + err_message)
	err_file.close()


def rename_reads(in_file_name):

	# This function renames reads in fastq files to make them distinguishable in direct and reverse readings

	with open(in_file_name, 'rU') as in_file:
		sys.stderr.write("Processing %s\n" % (in_file_name))
		out_file_name = in_file_name[:-6]+'.rn.fastq'
		out_file = open(out_file_name, 'w')

		for line in in_file:
			if line.startswith('@') and len(line.split(' ')) > 1:
				out_line = line.split(' ')[0] + line.split(' ')[1][0] + '\n'
			else:
				out_line = line
			out_file.write(out_line)
		out_file.close()
	return out_file_name




if __name__ == '__main__':
	args = parse_command_line_arguments()

	first_file = args.fastq_R1_file
	second_file = args.fastq_R2_file
	name = args.name
	name_ref = args.name_ref

	assert first_file.endswith('.fastq')
	assert second_file.endswith(".fastq")
	assert name_ref == 'canFam3' or 'bosTau7'

	first_rn_file = rename_reads(first_file) # rn = renamed
	second_rn_file = rename_reads(second_file)


	i = 0 # This variable will show number of command (0 - first, 1 - second and so on)

	command_list=	[
			['~/.local/bin/cutadapt -a AGATCGGAAGAGC -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3 -o - ' + first_rn_file],
			['~/.local/bin/cutadapt -a AGATCGGAAGAGC -a CCACATNNNNNNCTCGAGTCGG -g CCGACTCGAGNNNNNNATGTGG -n 3 -o - ' + second_rn_file],
			['~/programs/bowtie2-2.2.4/bowtie2 -p4 --local -x ' + name_ref + ' -1 ' + name + '.R1.ca.fastq' + ' -2 ' + name + '.R2.ca.fastq'],
			['~/programs/bowtie2-2.2.4/bowtie2 -p4 --local -x ' + 'hg19' + ' -1 ' + name + '.R1.ca.fastq' + ' -2 ' + name + '.R2.ca.fastq'],
			['samtools view -bS ' + name + ".ca." + name_ref + ".pe.sam"],
			['samtools sort ' + name + ".ca." + name_ref + ".pe.bam " + name + ".ca." + name_ref + ".pe.srt"],
			['samtools view -bS ' + name + ".ca.hg19.pe.sam"],
			['samtools sort ' + name + ".ca.hg19.pe.bam " + name + ".ca.hg19.pe.srt"],
			]

	while i < 8: # Through this cycle 7 commands will be executed
		creating_files(command_list[i], i)
		i += 1

	sys.stderr.write("Completed\n")
