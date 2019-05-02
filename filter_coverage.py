#!/usr/bin/env python
import sys
from Bio import SeqIO
import re, getopt

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"l:f:c:",[])
    if len(args) > 0:
        raise getopt.GetoptError("Please ensure that you have entered command correctly",None)
except getopt.GetoptError:
    print("Usage: python filter_coverage.py -l min_contig_length -f file_to_filter -c minimum_coverage")
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-l"):
        min_length = int(arg)
    elif opt in ("-f"):
        contigfile = str(arg)
    elif opt in ("-c"):
        min_coverage = float(arg)
    else:
        print("Unrecognized option: " + opt)


#cov_pattern = re.compile("cov_([0-9.]+)$")
cov_pattern = re.compile("cov_([0-9.]+)")


#min_coverage, fasta_file_path = sys.argv[1:]
#with open(contigfile.replace('fasta', 'filter{}cov.fasta'.format(min_coverage)), 'rUw') as filtered_fasta:
with open(contigfile.replace('.fasta', '_fullyfiltered.fasta'), 'w') as write_fasta:
    with open(contigfile.replace('.fasta', '_filtered.fasta'), 'rU') as filtered_fasta:
    


#with open(fasta_file_path, 'rU') as input_fasta:
        def filtered_contigs_generator(min):
            for contig in SeqIO.parse(filtered_fasta, 'fasta'):
                result = cov_pattern.search(contig.name)
                if result:
                    if float(result.group(1)) >= min:
                        #contig.name = contig.name.replace("NODE",contigfile.rstrip(".fasta").lstrip("./"))
                        yield contig
                
                
        SeqIO.write(filtered_contigs_generator(float(min_coverage)), write_fasta, 'fasta')
