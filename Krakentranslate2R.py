#!/usr/bin/env python

# Script that takes a krakentranslate file, and breaks down results sample-wise so as to allow easy plotting in R

#Need:
#- Separate file by sample
#- Separate contigs in sample by species/other phylogenetic level
#- Each contig must have a LEN and COV parameter

#Pseudocode:
#Open file
#Read line by line (generator)
#Enter each line into a dic with samples as keys
#Within each sample-dic, enter each species as keys
#Within each species-dic, Enter contigs like {1: {LEN: 25000, COV: 20.2}, 2 ...}


import sys
from os.path import basename

master_results = {}

if len(sys.argv) > 2:
	infile = sys.argv[1]
	outfile = sys.argv[2]
else:
	sys.exit("Usage: Krakentranslate2R.py infile.krakentranslate outfile.Rtable")

with open(infile, "r") as krakenfile:
	for line in krakenfile:
		# Split by tab into [0] = contig name [1] = annotation
		linesplit = line.split("\t")
		contigname = linesplit[0].split("_")
		lengthindex = contigname.index("length")
		covindex = contigname.index("cov")
		
		contignumber = contigname[lengthindex - 1]
		samplename = "_".join( [ contigname[x] for x in range(lengthindex-1) if contigname[x] != "NODE" ] )
		contiglength = contigname[lengthindex + 1]
		contigcov = contigname[covindex + 1]
		
		annotation = linesplit[1].split(";")
		
		if len(annotation) >= 9: # An identification at least to species level
			species = annotation[8].rstrip("\n")
		else:
			species = annotation[-1].rstrip("\n")
		
		if samplename in master_results:
			if species in master_results[samplename]:
				master_results[samplename][species][contignumber] = {"Len": contiglength, "Cov": contigcov}
			else:
				master_results[samplename][species] = {contignumber: {"Len": contiglength, "Cov": contigcov}}			
		else:
			master_results[samplename] = {species: {contignumber: {"Len": contiglength, "Cov": contigcov}}}

# Create a table/data frame to send to R
with open(outfile, "w") as rfile:
	rfile.write('"Sample","Species","Contig","Length","Coverage"\n')
	for sample in master_results:
		# Sort species so that the most common occur first, then the second most common etc
		species = list(master_results[sample].keys())
		species = sorted(species, key=species.count, reverse=True)
		for specie in species:
			contignumbers = list(master_results[sample][specie].keys())
			contignumbers = sorted( [int(x) for x in contignumbers] )
			for c in contignumbers:
				outline = [sample, specie, str(c), master_results[sample][specie][str(c)]["Len"], master_results[sample][specie][str(c)]["Cov"]]
				rfile.write(",".join('"' + item + '"' for item in outline) + "\n")
