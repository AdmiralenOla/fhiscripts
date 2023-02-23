#!/usr/bin/env python

'''
Script for reading in Kraken MPA classifications and generating aggregated sample outputs.
Usage: python generate_aggregated_krona_output.py <input1.mpa> <input2.mpa> ... <inputN.mpa>
'''

import csv, argparse, os
import pandas as pd

parser = argparse.ArgumentParser(description='Generate aggregated Krona output')
parser.add_argument("-i", "--input_files",nargs='+',help="Path to KRONA file(s). (Accepts multiple)", required=True)
parser.add_argument("-o", "--outfile", required=True,help="Path to output file.", default="Aggregated_Krona_output.tsv")
parser.add_argument("-l", "--level", help="Taxonomic level.",required=True,type=str,choices=['Domain','Phylum','Class','Order','Family','Genus','Species'])
args = parser.parse_args()

def find_level(txstring,level):
	txsplit = txstring.split("|")
	l = txsplit[-1]
	if level == "Domain" and l.startswith("d__"):
		return l[3:]
	elif level == "Phylum" and l.startswith("p__"):
		return l[3:]
	elif level == "Class" and l.startswith("c__"):
		return l[3:]
	elif level == "Order" and l.startswith("o__"):
		return l[3:]
	elif level == "Family" and l.startswith("f__"):
		return l[3:]
	elif level == "Genus" and l.startswith("g__"):
		return l[3:]
	elif level == "Species" and l.startswith("s__"):
		return l[3:]
	else:
		return None

def main():
	results = {}
	for infile in args.input_files:
		# Call from data dir so that first part of path is samplename
		samplename = infile.split("/")[0]
		with open (infile,"r") as my_infile:
			text = csv.reader(my_infile,delimiter="\t")
			sampleresults = {}
			for line in text:
				taxa_string = line[0]
				taxa = find_level(taxa_string,args.level)
				count = line[-1]
				if taxa is not None:
					sampleresults[taxa] = count
			print("Number of counts in file %s: %s" % (samplename, str(len(sampleresults))))
			# NOTE: This will truncate so that only index keys are used
			results[samplename] = pd.Series(sampleresults)
	resultsdf = pd.DataFrame(results).fillna(0)
	with open(args.outfile,'w') as outfile:
		resultsdf.to_csv(outfile, sep="\t")

if __name__ == '__main__':
    main()