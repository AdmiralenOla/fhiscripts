#!/bin/bash
# Script to run trimmomatic and spades on a run
# Must be in Run folder. (All isolates exist as folders in the current directory)

# USAGE: Trim_and_assemble.sh

VERSION="1.0"

basedir=$(pwd)
runname=${basedir##*/}

for dir in $(ls -d */)
do
	cd ${dir}
	strain=${PWD##*/}
	# Find R1 and R2, get results into an array
	R1=($(find -name "*R1*"))
	# If more than 1 R1 in dir, (NextSeq reads), merge reads together
	if [ ${#R1[@]} -gt 1 ]; then		
		Rlen=${#R1[@]}
		for (( i=0; i<${Rlen}; i++ ));
		do
			cat ${R1[$i]} >> ${strain}_merged_R1.fastq.gz
			cat ${R1[$i]/R1/R2} >> ${strain}_merged_R2.fastq.gz
		done
		R1=($(find -name "*merged_R1.fastq.gz"))
	fi

	#R1=$(ls *R1*)
	#NOTE - If trimmed results already exist - Don't do trimmomatic
	check=($(find -name "*1P.fastq.gz"))
	if ! [ ${#check[@]} -gt 0 ]; then
		trimmomatic PE -basein ${R1[0]} -baseout ${R1[0]%%_R1*}.fastq.gz ILLUMINACLIP:/opt/conda/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36
	fi

	newR1=($(find -name "*1P.fastq.gz"))
	newR2=($(find -name "*2P.fastq.gz"))

	if ! [ -f Output ] || [ -f ${basedir}/Spades_assembly/${strain} ]; then
		spades.py -o Output --careful --cov-cutoff auto -t 4 -1 ${newR1[0]} -2 ${newR2[0]}
	fi

	cd "${basedir}"
done

cd ${basedir}

# Copy all contig files to new directory and name appropriately
spades_output_dir="./Spades_assembly"
if ! test -f "Spades_assembly"; then
	mkdir "./Spades_assembly"	
fi

for dir in $(ls -d */)
do
	#cp ${dir}Data/Intensities/BaseCalls/contigs.fasta ./Spades_assembly/${dir%-????????/}.fasta
	# IF FILES DONT HAVE THE SILLY -????? etc suffix:
	if test -f ${dir}Output/contigs.fasta; then
		cp ${dir}Output/contigs.fasta "${spades_output_dir}/${dir%/}.fasta"
	fi
done

# Run filtration
echo ""
echo "Filtering bad contigs"
echo ""
cd "${spades_output_dir}"
filter_bad_contigs.sh -d ${spades_output_dir}

# Run Kraken on unfiltered to screen for contamination
#cd Filtered
mkdir Kraken
KRAKENDB="/opt/minikraken/minikraken"
#export KRAKEN_DEFAULT_DB="/home/ngs1/miniconda3/share/kraken/minikraken_20171019_8GB" # CHANGE DB
#export KRAKEN_NUM_THREADS=4
#kraken --fasta-input --preload --threads 4 --output Kraken/All_contigs.kraken *.fasta #UPDATE FOR KRAKEN2
#kraken-report Kraken/All_contigs.kraken > Kraken/All_contigs.krakenreport 
#kraken-translate Kraken/All_contigs.kraken > Kraken/All_contigs.krakentranslate
#rm Kraken/All_contigs.kraken

# Re-run Kraken on filtered to make sure contigs are good
echo ""
echo "Running Kraken"
echo ""
cd Filtered
mkdir Kraken
kraken2 --db "$KRAKENDB" --threads 4 --report Kraken/All_contigs.krakenreport --use-names --output Kraken/All_contigs.kraken *.fasta # UPDATE FOR KRAKEN2
#kraken-report Kraken/All_contigs.kraken > Kraken/All_contigs.krakenreport
#kraken-translate Kraken/All_contigs.kraken > Kraken/All_contigs.krakentranslate
rm Kraken/All_contigs.kraken

# Find PhiX contigs and exclude
# < <( does process substitution
echo ""
echo "Removing phiX contigs"
echo ""

readarray phiX < <(grep "phiX" Kraken/All_contigs.krakentranslate)

# Create method to find substring index
# Returns -1 if not found. Returns index (0-based) if substring b found in a
# Call like
# strindex "$a" "$b"
# Warning: If multiple occurences, not guaranteed to find first/last match
strindex() {
	x="${1%%$2*}"
	[[ $x = $1 ]] && echo -1 || echo ${#x}
}

for elem in "${phiX[@]}"
do
	separatorindex=$(strindex "$elem" "_length")
	# This index separates the Name and contig number from the rest of the contig name, such as _length_ etc
	nameandcontig="${elem:0:$separatorindex}"
	# Now we have a string with the format NameofStrain_Contignumber
	
	contignameindex=$(strindex "$elem" "root")
	contigname="${elem:0:$((contignameindex - 1))}"
	
	# Open FASTA file and remove contig
	filename="${nameandcontig%_*}.fasta"
	
	# Need to export contigname so that perl can see it in its $ENV variable
	export contigname
	perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/} $ENV{contigname} }print if not $c' $filename > ${filename}_nophix
	mv ${filename}_nophix ${filename}
done

# Run Krakentranslate2R for Affyboy-style figures

Krakentranslate2R.py Kraken/All_contigs.krakentranslate Kraken/All_contigs.Rtable

Rscript --vanilla /usr/bin/fhiscripts/Rscript_Kraken.R "$basedir"/"${spades_output_dir}"/Filtered/Kraken/ "$basedir"/"${spades_output_dir}"/Filtered/Kraken/All_contigs.Rtable
	
# REMOVE ALL UNNECESSARY FILES

cd "$basedir"

echo "You are now in $(pwd). Remove all temporary non-result files (e.g. spades intermediary files)?"
echo "[y]es / [n]o"

read yesnodel


if [ "$yesnodel" == "y" ]; then
	rm -rf */Output
	rm */*1P.fastq.gz
	rm */*1U.fastq.gz
	rm */*2P.fastq.gz
	rm */*2U.fastq.gz
else
	echo "Exiting without removing intermediary files"
	exit
fi