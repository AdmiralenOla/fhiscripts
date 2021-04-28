#!/usr/bin/env python

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re, argparse

def read_FASTA_Bio(f):
	"""Read a FASTA file in as a SeqIO"""
	return SeqIO.read(f, "fasta")

def read_MUM_SNP(f):
	'''Read a Mummer SNP file'''
	# Name: x1.snps.txt
	change = {} # Structure : { 1 : {Refpos: 19, Refbase: "A", Newbase: "C"}} 
	i = 1
	pattern = "^\s+(\d+)\s+(A|C|G|T|.) (A|C|G|T|.)"
	with open(f, 'rU') as mumfile:
		for line in mumfile:
			checkmatch = re.match(pattern, line)
			if checkmatch:
				if checkmatch[2] != "." and checkmatch[3] != ".":
					change[i] = {"Refpos": checkmatch[1], "Refbase": checkmatch[2], "Newbase": checkmatch[3]}
					i += 1
	return change

def main():
	'''
	Takes in a reference chromosome (must be complete in a single contig) and a list of mummer snps.txt files. Creates an alignment where the SNPs are patched into each strain
	'''
	parser = argparse.ArgumentParser(description='Takes in a reference chromosome (must be complete in a single contig) and a list of mummer snps.txt files. Creates an alignment where the SNPs are patched into each strain')
	parser.add_argument('--reference', help='A single-contig FASTA file, serving as the reference. Contig header must be name of reference')
	parser.add_argument("--outfile", help='The name of the output file (Default: out.aln)', default='out.aln')
	parser.add_argument('mumsnp', nargs='+', metavar='MUM_SNP', help='MUMMER SNP file(s) that should be added to the alignment. Multiple allowed. Must be of style filename.snps.txt')
	args = parser.parse_args()

	Ref = read_FASTA_Bio(args.reference)
	Alignment = {1:Ref}
	for M in args.mumsnp:
		# Open SNP file
		# Create copy of reference
		# For each SNP substitution (NOT INDELS) - replace reference base with isolate base
		# Patch into multialignment (AlignIO ?)
		name = re.search("(\w+)\.snps\.txt", M).group(1)
		changes = read_MUM_SNP(M)
		
		newseq = Seq(''.join([b for b in Ref.seq]), IUPAC.unambiguous_dna).tomutable()
		for change in changes:
			if Ref.seq[int(changes[change]["Refpos"])-1] == changes[change]["Refbase"]:
				newseq[int(changes[change]["Refpos"])-1] = changes[change]["Newbase"]
			else:
				print("Something wrong")
		newseq = newseq.toseq()
		Iso = SeqRecord(id=name,description="",name="",seq=newseq)
		Alignment[len(Alignment) + 1] = Iso
	SeqIO.write([v for v in Alignment.values()], args.outfile, "fasta")

if __name__ == '__main__':
	main()
