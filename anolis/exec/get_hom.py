#!/usr/bin/env python


from collections import namedtuple, defaultdict, OrderedDict
import sys
import argparse


def parse_command_line_arguments():

	parser = argparse.ArgumentParser(description=	
					"""Outputs homologies and their scores for each scaffold
					from .net file
					"""
					)
	parser.add_argument("net_file", help="net file (.net)")

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_command_line_arguments()

	top_net = defaultdict(list)
	with open(args.net_file) as input_file:
		for line in input_file:
			LL = line.split()
			if LL[0] == 'net': # new chromosome
				t_chr = LL[1]
			elif LL[16] == 'top': # top level of net
				top_net[t_chr].append( [LL[3],LL[10]] ) # Q_chr, alignment score
	for (t_chr,fill) in top_net.items():
		aln_scores = defaultdict(int)
		for F in fill:
			aln_scores[F[0]] += int(F[1]) # sum up scores for Q_chr
		SAS = OrderedDict(sorted(aln_scores.items(), key=lambda t: t[1], reverse=True)) # sort by score
		SAL = []
		for L in SAS.items():
			SAL+=list(L)
		SAL = [str(L) for L in SAL]
		print t_chr + '\t' +','.join(SAL)

	sys.stderr.write("Complete!\n")
