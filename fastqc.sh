#!/bin/bash

# Script for calculating the coverage of an assembly. Requires input: FASTA file (contigs), FASTQ F/R
# Uses coverM (https://github.com/wwood/coverm)
eval "$(conda shell.bash hook)"
conda activate qc

# DEFAULT ARGUMENTS

FWD=false
REV=false

programname=$0

function usage {
	echo "usage $programname -f FWD.fastq -r REV.fastq"
	echo "  -f		Path to forward (1) reads in FASTQ (.gz) format"
	echo "  -r		Path to reverse (2) reads in FASTQ (.gz) format"
	exit 1
}

while getopts f:r: option
do
	case "${option}"
		in
		f) FWD=${OPTARG};;
		r) REV=${OPTARG};;
	esac
done

# VERIFY INPUT


if ! test -f $FWD; then
	echo "Forward FASTQ file not found. Exiting"
	usage
fi

if ! test -f $REV; then
	echo "Reverse FASTQ not found. Exiting"
	usage
fi

fastqc $FWD $REV

conda deactivate