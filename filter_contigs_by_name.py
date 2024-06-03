#!/usr/bin/env python
# Usage: python filter_contigs.py min_length contig_file_path
import sys, argparse, csv
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description='Filter contigs / alignment by name')
    parser.add_argument("-i", "--input_file",help="Path to input file in FASTA format", required=True)
    parser.add_argument("-f", "--filter_file", help="Path to file containing contigs to be included. One ID per line")
    parser.add_argument("-o", "--out_file", help="Path to output file.", required=True)
    args = parser.parse_args()

    with open(args.filter_file,'r') as filterfile:
        filter_list = [x.rstrip("\n") for x in filterfile.readlines()]

    with open(args.out_file,'w') as outfile:
        with open(args.input_file,'r') as infile:
            def filtered_contigs_generator():
                for contig in SeqIO.parse(infile,'fasta'):
                    if contig.id in filter_list:
                        yield contig
            SeqIO.write(filtered_contigs_generator(), outfile, 'fasta')


if __name__ == "__main__":
    main()