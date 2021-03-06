#!/bin/bash
# Script to run trimmomatic and spades on a run
# Must be in Run folder. (All isolates exist as folders in the current directory)

# USAGE: Trim_and_assemble.sh

VERSION="1.3"

basedir=$(pwd)
runname=${basedir##*/}
KRAKENDB="/opt/minikraken/minikraken"

# Version 1.2 - Array to hold overview of files for each strain
declare -A my_array

for dir in $(ls -d */)
do
	# Skip dir if name equal to Spades_assembly
	[ $dir = "Spades_assembly/" ] && continue
	# Skip dir if name equal to QC
	[ $dir = "QC/" ] && continue
	# Enter dir
	cd ${dir}
	# Get strain name
	strain=${PWD##*/}
	# Find R1 and R2, get results into an array
	R1=($(find -name "*R1*fastq.gz" | sort ))
	# If more than 1 R1 in dir, (NextSeq reads), merge reads together
	# However, if merged reads already exists, skip this step
	R1merged=($(find -name "*merged*fastq.gz"))
	if [ ${#R1[@]} -gt 1 ] && [ ${#R1merged[@]} -eq 0 ]; then		
		Rlen=${#R1[@]}
		echo "Merging ${Rlen} read files in ${strain}"
		for (( i=0; i<${Rlen}; i++ ));
		do
			cat ${R1[$i]} >> ${strain}_merged_R1_001.fastq.gz
			cat ${R1[$i]/R1/R2} >> ${strain}_merged_R2_001.fastq.gz
		done
		# Set R1 to be the merged file
		R1=($(find -name "*merged_R1_001.fastq.gz"))
	fi

	#R1=$(ls *R1*)
	#NOTE - If trimmed results already exist - Don't do trimmomatic
	check=($(find -name "*1P.fastq.gz"))
	if ! [ ${#check[@]} -gt 0 ]; then
		echo "Trimming strain ${strain}"
		trimmomatic PE -basein ${R1[0]} -baseout ${R1[0]%%_R1_001.fastq.gz}.fastq.gz ILLUMINACLIP:/opt/conda/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:15 MINLEN:36
	fi

	# newR1/R2 should be the trimmed files - Note these paths are relative to the 
	newR1=($(find -name "*1P.fastq.gz"))
	newR2=($(find -name "*2P.fastq.gz"))

	# Add to my_array
	my_array["${strain}_R1"]="${basedir}/${dir}${newR1}"
	my_array["${strain}_R2"]="${basedir}/${dir}${newR2}"

	# Run spades unless fasta file already exists in Spades_assembly OR dir Output exists in current folder
	echo "Running Spades on strain ${strain}"
	if ! ( [ -f Output ] || [ -f ${basedir}/Spades_assembly/${strain}.fasta ] ); then
		spades.py -o Output --isolate --cov-cutoff auto -t 4 -1 ${newR1[0]} -2 ${newR2[0]}
	fi

	cd "${basedir}"
done

cd ${basedir}

# Copy all contig files to new directory and name appropriately
spades_output_dir="${basedir}/Spades_assembly"
if ! test -d "Spades_assembly"; then
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
filter_bad_contigs.sh

# Re-run Kraken on filtered to make sure contigs are good
echo ""
echo "Running Kraken on unfiltered contigs"
echo ""

# Run Kraken on unfiltered to screen for contamination
#cd Filtered
if ! test -d "Kraken"; then
	mkdir Kraken
fi
#export KRAKEN_DEFAULT_DB="/home/ngs1/miniconda3/share/kraken/minikraken_20171019_8GB" # CHANGE DB
#export KRAKEN_NUM_THREADS=4
kraken --db "$KRAKENDB" --fasta-input --preload --threads 4 --output Kraken/All_contigs.kraken *.fasta #UPDATE FOR KRAKEN2
kraken-report --db "$KRAKENDB" Kraken/All_contigs.kraken > Kraken/All_contigs.krakenreport 
kraken-translate --db "$KRAKENDB" Kraken/All_contigs.kraken > Kraken/All_contigs.krakentranslate
rm Kraken/All_contigs.kraken

# Find PhiX contigs and exclude
# Note that this needs to operate on UNFILTERED contigs

readarray phiX < <(grep "phiX" Kraken/All_contigs.krakentranslate)
# Can change into filtered dir after finding which contigs are phiX
cd Filtered


# < <( does process substitution
echo ""
echo "Removing phiX contigs from filtered contigs"
echo ""

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

# Re-run Kraken on filtered to make sure contigs are good
echo ""
echo "Running Kraken on filtered contigs"
echo ""
if ! test -d "Kraken"; then
	mkdir "Kraken"
fi
kraken --db "$KRAKENDB" --threads 4 --fasta-input --output Kraken/All_contigs.kraken --preload *.fasta
kraken-report --db "$KRAKENDB" Kraken/All_contigs.kraken > Kraken/All_contigs.krakenreport
kraken-translate --db "$KRAKENDB" Kraken/All_contigs.kraken > Kraken/All_contigs.krakentranslate
rm Kraken/All_contigs.kraken


# Run Krakentranslate2R for Affyboy-style figures

Krakentranslate2R.py Kraken/All_contigs.krakentranslate Kraken/All_contigs.Rtable

Rscript --vanilla /usr/bin/Rscript_Kraken.R "${spades_output_dir}"/Filtered/Kraken/ "${spades_output_dir}"/Filtered/Kraken/All_contigs.Rtable

# Version 1.2 - DIAGNOSTICS PER STRAIN
echo "Creating run diagnostics per genome"
echo "Strain Ncontigs Length Coverage" >> "${spades_output_dir}/Spades_FHI_RUNDIAGNOSTICS.txt"
for output in *.fasta
do
	strainname=${output%.fasta}
	this_R1="${my_array[${strainname}_R1]}"
	this_R2="${my_array[${strainname}_R2]}"
	echo "Calculating coverage of strain ${strainname}, fastafile ${output}, R1 ${this_R1}, R2 ${this_R2}"
	echo "${strainname} $(coverage_calculator.sh -c ${output} -f ${this_R1} -r ${this_R2})" >> "${spades_output_dir}/Spades_FHI_RUNDIAGNOSTICS.txt"
done


# TAG WITH VERSION AND DATE
echo "VERSION=${VERSION}" >> "${spades_output_dir}/Spades_FHI_RUNINFO.txt"
echo "DATE=$(date)" >> "${spades_output_dir}/Spades_FHI_RUNINFO.txt"


# REMOVE ALL UNNECESSARY FILES


echo "You are now in $(pwd). Remove all temporary non-result files (e.g. spades intermediary files)?"
echo "[y]es / [n]o"

read yesnodel


if [ "$yesnodel" == "y" ]; then
	rm -rf */Output
	rm */*1P.fastq.gz
	rm */*1U.fastq.gz
	rm */*2P.fastq.gz
	rm */*2U.fastq.gz
	rm */*.kraken
else
	echo "Exiting without removing intermediary files"
	exit
fi
