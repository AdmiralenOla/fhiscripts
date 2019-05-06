#!/usr/bin/env python

'''Rename contigs of a FASTA file to replace NODE with FILENAME (sans .fasta) '''

import argparse, os


def main():
    '''Execute renaming.'''

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Rename FASTA files.', epilog='Work out those contigs.')
    parser.add_argument('-i', '--input', help='indicate input FASTA file', required=True)
    parser.add_argument('--pre', help='string pre contig count', type=str, default='')
    parser.add_argument('--pos', help='string post contig count', type=str, default='')
    parser.add_argument('-o', '--output', help='indicate output FASTA file', required=False)
    parser.add_argument('-f','--full', dest="full", help='should the contig name be just the isolate name (BIGSdb upload)', action="store_true")
    parser.add_argument('-nf', '-notfull', dest="full", action="store_false")
    parser.set_defaults(full=True)
    args = parser.parse_args()

    # Open FASTA.
    fasta_in = open(args.input)
    strain_name = os.path.splitext(args.input)[0]
    strain_name = strain_name.lstrip("./Filtered/")

    # Create FASTA output file.
    if args.output:
        rename = False
        fasta_out = open(args.output, 'w')
    else:
        rename = True
        fasta_out = open(args.input + "_out", "w")

    # Start counter.
    #count = 1
    
    # Parse file and write to output.
    print('Parsing %s...' % args.input)
    for line in fasta_in.readlines():
        if line.startswith('>'):
            firstunderscore = line.find("_")
            
            if bool(args.full):
                contig_id = '>' + strain_name + '\n'
            else:
                contig_id = '>' + args.pre + strain_name + line[firstunderscore:]
                
            fasta_out.write(contig_id)
            #count += 1
        else:
            fasta_out.write(line)

    # Finish.
    fasta_out.close()
    fasta_in.close()
    #print('Wrote %d contigs to %s.' % (count, args.output))
    
    if rename:
        os.rename(args.input + "_out", args.input)

if __name__ == '__main__':
    main()
