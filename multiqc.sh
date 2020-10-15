#!/bin/bash

# Script for calculating the coverage of an assembly. Requires input: FASTA file (contigs), FASTQ F/R
# Uses coverM (https://github.com/wwood/coverm)
eval "$(conda shell.bash hook)"
conda activate qc

multiqc .

conda deactivate
