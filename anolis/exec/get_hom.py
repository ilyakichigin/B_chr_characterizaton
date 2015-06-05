#!/usr/bin/env python


from collections import namedtuple, defaultdict, Counter, OrderedDict
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

	File = args.net_file

	TopNet = defaultdict(list)
	with open(File) as f:
		for Line in f:
			LL = Line.split()
			if LL[0] == 'net': # new chromosome
				T_chr = LL[1]
			elif LL[16] == 'top': # top level of net
				TopNet[T_chr].append( [LL[3],LL[10]] ) # Q_chr, alignment score
	for (T_chr,Fill) in TopNet.items():
		AlnScores = defaultdict(int)
		for F in Fill:
			AlnScores[F[0]] += int(F[1]) # sum up scores for Q_chr
		SAS = OrderedDict(sorted(AlnScores.items(), key=lambda t: t[1], reverse=True)) # sort by score
		SAL = []
		for L in SAS.items():
			SAL+=list(L)
		SAL = [str(L) for L in SAL]
		print T_chr + '\t' +','.join(SAL)
