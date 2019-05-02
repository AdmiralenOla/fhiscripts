#!/bin/sh

# This script automizes the filtration of multi-contig fasta files to exclude contigs with length <= L and coverage <= C

# Default values:

# USAGE: filter_bad_contigs.sh -l minimumcontiglength -c minimumcoverage -d directoryofcontigfiles

echo "Filtering contigs in directory"
LENGTH=500
COVERAGE=2.0
DIRECTORY="./"
FULLWIPE=false

for i in "$@"
do
	
case $i in
	-l=*|--length=*)
	LENGTH="${i#*=}"
	shift # past-argument = value
	;;
	-c=*|--coverage=*)
	COVERAGE="${i#*=}"
	shift # past-argument = value
	;;
	-d=*|--directory=*)
	DIRECTORY="${i#*=}"
	shift
	;;
	-f)
	FULLWIPE=true
	shift
	;;
	*)
	# unknown option
	echo "Unknown argument ${i}"
	;;
esac
done

echo "LENGTH = ${LENGTH}"
echo "COVERAGE = ${COVERAGE}"
echo "DIRECTORY = ${DIRECTORY}"
if [ "$FULLWIPE" = true ]; then
	echo "WIPING CONTIG NAMES AND PREPARING FOR BIGSDB"
fi


# SEND FILE NOT DIRECTORY
all=$(ls ${DIRECTORY}*.fasta)
mkdir $DIRECTORY\Filtered
for x in $all
do
	filter_contigs.py -l $LENGTH -f $x -c $COVERAGE
	filter_coverage.py -l $LENGTH -f $x -c $COVERAGE
	#mv ${x%.fasta}_fullyfiltered.fasta ${x%.fasta}_filtered.fasta
	#mv ${x%.fasta}_filtered.fasta $DIRECTORY\Filtered/
	mv ${x%.fasta}_fullyfiltered.fasta $DIRECTORY\Filtered/$x
	rm ${x%.fasta}_filtered.fasta
	if [ "$FULLWIPE" = true ]; then
		rename_contigs.py -i $DIRECTORY\Filtered/${x#./} -f
	else
		rename_contigs.py -i ${DIRECTORY}${x#./} -nf
		rename_contigs.py -i $DIRECTORY\Filtered/${x#./} -nf
	fi
done
