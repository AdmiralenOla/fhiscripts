#!/usr/bin/env python
# Usage: python filter_contigs.py min_length contig_file_path
import sys, getopt
from Bio import SeqIO

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"l:f:c:",[])
    if len(args) > 0:
        raise getopt.GetoptError("Please ensure that you have entered command correctly",None)
except getopt.GetoptError:
    print("Usage: python filter_contigs.py -l min_contig_length -f file_to_filter -c minimum_coverage")
    sys.exit(2)
    
for opt, arg in opts:
    if opt in ("-l"):
        min_length = int(arg)
    elif opt in ("-f"):
        contigfile = str(arg)
    elif opt in ("-c"):
        mincoverage = arg
    else:
        print("Unrecognized option: " + opt)


#min_length, fasta_file_path = sys.argv[1:]
#with open(contigfile.replace('fasta', 'filter{}.fasta'.format(min_length)), 'w') as filtered_fasta:
with open(contigfile.replace('.fasta', '_filtered.fasta'), 'w') as filtered_fasta:
    with open(contigfile, 'rU') as input_fasta:
        def filtered_contigs_generator(min):
            for contig in SeqIO.parse(input_fasta, 'fasta'):
                if len(contig) >= min:
                    yield contig
        
        SeqIO.write(filtered_contigs_generator(int(min_length)), filtered_fasta, 'fasta')
        #SeqIO.write(filtered_contigs_generator(int(min_length)), sys.stdout, 'fasta')
