#!/bin/bash

# Script for calculating the coverage of an assembly. Requires input: FASTA file (contigs), FASTQ F/R
# Uses coverM (https://github.com/wwood/coverm)
eval "$(conda shell.bash hook)"
conda activate coverage

# DEFAULT ARGUMENTS

REF=false
FWD=false
REV=false

programname=$0

function usage {
	echo "usage $programname -c contigs.fasta -f FWD.fastq -r REV.fastq"
	echo "  -c		Path to reference FASTA file (contigs)"
	echo "  -f		Path to forward (1) reads in FASTQ (.gz) format"
	echo "  -r		Path to reverse (2) reads in FASTQ (.gz) format"
	exit 1
}

function stats {
	awk '/^>/ { seqtotal+=seqlen
            seqlen=0
            seq+=1
            next
            }
    {
    seqlen += length($0)
    }     
    END{
        print "Contigs: " seq " Length: " seqtotal+seqlen
    }' $1
}

while getopts c:f:r: option
do
	case "${option}"
		in
		c) REF=${OPTARG};;
		f) FWD=${OPTARG};;
		r) REV=${OPTARG};;
	esac
done

# VERIFY INPUT

if ! test -f $REF; then
	echo "Contigs file not found. Exiting"
	usage
fi

if ! test -f $FWD; then
	echo "Forward FASTQ file not found. Exiting"
	usage
fi

if ! test -f $REV; then
	echo "Reverse FASTQ not found. Exiting"
	usage
fi

COV=$(coverm genome -1 $FWD -2 $REV -r $REF -t 4 --single-genome -m mean -q)

STATS=$(stats ${REF})

echo $STATS $COV | cut -d " " -f 2,4,9

conda deactivate