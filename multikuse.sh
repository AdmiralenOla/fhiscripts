#!/bin/bash
# Script to run trimmomatic and spades on a run
# Must be in Run folder. (All isolates exist as folders in the current directory)

# USAGE: Trim_and_assemble.sh

VERSION="1.0"

basedir=$(pwd)
runname=${basedir##*/}
KRAKENDB="/opt/minikraken/minikraken"

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
	R1=($(find -name "*R1*" | sort ))
	# If more than 1 R1 in dir, (NextSeq reads), merge reads together
	# However, if merged reads already exists, skip this step
	R1merged=($(find -name "*merged*"))
	if [ ${#R1[@]} -gt 1 ] && [ ${#R1merged[@]} -eq 0 ]; then		
		# Set R1 to be the merged file
		R1=($(find -name "*merged_R1_001.fastq.gz"))
	fi

	#R1=$(ls *R1*)
	#NOTE - If trimmed results already exist - Don't do trimmomatic
	check=($(find -name "*1P.fastq.gz"))
	if [ ${#check[@]} -gt 0 ]; then
		R1=($(find -name "*1P.fastq.gz"))
		R2=($(find -name "*2P.fastq.gz"))
	else
		R2="${R1/R1/R2}"
	fi

	# Run FASTQC
	fastqc.sh -f $R1 -r $R2

	# Run kraken on samples
	echo "Running Kraken on strain ${strain}"
	if ! test -f ${strain}.kraken; then
		kraken --db "$KRAKENDB" --threads 4 --fastq-input --gzip-compressed --paired --output ${strain}.kraken --preload $R1 $R2
	fi
	# Kraken-report
	if ! test -f ${strain}.strainreport; then
		kraken-report --db "$KRAKENDB" ${strain}.kraken > ${strain}.strainreport
	fi

	cd "${basedir}"
done

cd ${basedir}

# Run multiqc
cd "$basedir"
all_strains=$(ls -d */)
delete=("Spades_assembly/" "QC/")
for del in ${delete[@]}
do
	all_strains="${all_strains[@]/$del}"
done
multiqc.sh 
