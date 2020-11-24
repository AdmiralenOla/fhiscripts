#!/usr/bin/env python

'''
mask_blast_hit_with_N.py

Author: Ola Brynildsrud
Date Started: 20/11/2020

Accepts a query sequence in FASTA format and a subject sequence in
FASTA format. The coordinates of the query in the subject is masked with N
'''
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Alphabet import IUPAC, SingleLetterAlphabet
from Bio.Blast import NCBIXML
from pkg_resources import resource_string, resource_filename
#from .__init__ import __version__
import argparse, time, os, sys
__version__ = '0.2b'

def rev_comp(seq):
    """Returns the reverse complement of a given sequence"""
    function = {"A":"T", "T":"A", "C":"G", "G":"C"}
    output = [function[x] for x in seq[::-1]]
    return "".join(output)

def rev_comp_Bio(seq):
    """Takes SeqIO object and reverse complements"""
    return seq.reverse_complement(id=seq.id, description= "reverse complemented")

def rotate_seq(seq,coord):
    """Rotates the circular sequence so that coord is number 1. Should be able to handle Seq objects"""
    #Make sure coord is 0-based
    coord -= 1
    return seq[coord:] + seq[:coord]

def mask_seq(seq, start, end, id, length):
    '''Replaces the coords between start and end with N'''
    masked = []
    insert = Seq("N" * length, SingleLetterAlphabet())
    # Find contig
    for seq_record in seq:
        if seq_record.id == id:
            newseq = seq_record.seq[:start-1] + insert + seq_record.seq[end:]
            #seq_record.seq = newseq
            masked.append(SeqRecord(id=seq_record.id, description=seq_record.description,seq=newseq))
        else:
            masked.append(seq_record)
    return masked
    # Mask coords

def read_FASTA(f):
    """Read in a file in FASTA format"""
    print ("File: ", f)
    print ("Reading in FASTA...", end=' ')

    in_file = open(f,'r')
    seqDict = {}
    name = None
    for line in in_file:
        line = line.rstrip()
        if line[0]=='>':
            name = line[1:]
            seqDict[name] = ''
        else:
            seqDict[name] = seqDict[name] + line

    print ("DONE!")
    return seqDict

def find_island(seq, db, minlength, minperc):
    """Find the location of GGI and return the coordinates. Allow multiple hits"""
    blastn_cline = NcbiblastnCommandline(query=seq, subject=db, outfmt=5, out='GGI.xml')
    print("Running BLAST...")
    blastn_cline()
    blast_record = NCBIXML.parse(open('GGI.xml'))
    hits = []
    nohit = {'title': 'None', 'score': 0.0, 'identities': 0, 'frame':(0,0), 'gaps':0, 'positives':0, 'length':0, 'subject_start':0, 'subject_end': 0,'frame':(0,0)}
    hits.append(nohit)
    #highscore = 0.0
    # Min-length = 
    for Blast in blast_record:
        for al in Blast.alignments:
            for hsp in al.hsps:
                if hsp.align_length >= minlength:
                    percid = hsp.identities / hsp.align_length
                    if percid >= minperc:
                        hit = {'title': al.hit_id, 'score': hsp.score, 'identities': hsp.identities, 'frame':hsp.frame, 'gaps':hsp.gaps, 'positives':hsp.positives, 'length':hsp.align_length, 'subject_start':hsp.sbjct_start,'subject_end': hsp.sbjct_end,'frame':hsp.frame}
                        hits.append(hit)
                    
    return hits

def read_FASTA_Bio(f):
    """Read a FASTA file in as a SeqIO"""
    return SeqIO.parse(open(f,'rU'), "fasta")

def writeFile(fasta, outfile):
    """Writes the rotated FASTA to file"""
    print("Writing to file %s" % outfile)
    with open(outfile,'w') as f:
        for rec in fasta:
            try:
                SeqIO.write(rec,f,"fasta")
            except TypeError:
                print("ERROR: %s " % rec)

def main():
    """Main function of originator"""
    start = time.clock()
    parser = argparse.ArgumentParser(description='mask_blast_hit_with_N')
    parser.add_argument("-v", "--version", help="Installed version", action="version",
                        version="%(prog)s " + str(__version__))
    parser.add_argument("--query", help="Query FASTA file to BLAST for")
    parser.add_argument("fastaFile", help="Subject FASTA file to BLAST against")
    parser.add_argument("--outfile", help="Name of file to write masked FASTA to")
    parser.add_argument("--minlength", help="Minimum length of hit to mask. Default=100",default=100)
    parser.add_argument("--minidentity", help="Minimum percent identity to query to mask, over length of hit. Default=0.9",default=0.9)
    args = parser.parse_args()
    if args.outfile is None:
        args.outfile = args.fastaFile + ".masked.fasta"

    hits = find_island(args.query, args.fastaFile, args.minlength, args.minidentity)

    input_fasta = read_FASTA_Bio(args.fastaFile)
    #input_fasta = read_FASTA(args.fastaFile)
    if len(hits) == 1 and hits['title'][0] == "None":
        print("Unable to find any trace of the GGI island")
        writeFile(input_fasta,args.outfile)
    else:
        # MASK SEQUENCE
        masked_fasta = input_fasta
        for hit in hits[1:]:
            if hit['frame'] == (1,1):
                masked_fasta = mask_seq(masked_fasta, hit['subject_start'], hit['subject_end'], hit['title'], hit['length'])
            elif hit['frame'] == (1,-1):
                masked_fasta = mask_seq(masked_fasta, hit['subject_end'], hit['subject_start'], hit['title'], hit['length'])
            else:
                print("Something wrong with frame. Contact Ola")
                sys.exit(-1)
        writeFile(masked_fasta, args.outfile)


if __name__ == '__main__':
    main()
