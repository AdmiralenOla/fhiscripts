#!/usr/bin/env python

'''
Script for reading in Krona classifications and generating aggregated sample outputs.
Usage: python generate_aggregated_krona_output.py <input1.krona> <input2.krona> ... <inputN.krona>
'''

import csv, argparse, os
import pandas as pd

parser = argparse.ArgumentParser(description='Generate aggregated Krona output')
parser.add_argument("-i", "--input_files",nargs='+',help="Path to KRONA file(s). (Accepts multiple)", required=True)
parser.add_argument("-o", "--outfile", help="Path to output file.", default="Aggregated_Krona_output.tsv")
args = parser.parse_args()

def main():
	#results = {}
	results = {}
	for infile in args.input_files:
		# Call from data dir so that first part of path is samplename
		samplename = infile.split("/")[0]
		with open (infile,"r") as my_infile:
			text = csv.reader(my_infile,delimiter="\t")
			sampleresults = {}
			for line in text:
				taxa = line[-1]
				count = line[0]
				sampleresults[taxa] = count
			print(len(sampleresults))
			# NOTE: This will truncate so that only index keys are used
			results[samplename] = pd.Series(sampleresults)
	resultsdf = pd.DataFrame(results).fillna(0)
	with open(args.outfile,'w') as outfile:
		resultsdf.to_csv(outfile, sep="\t")

if __name__ == '__main__':
    main()