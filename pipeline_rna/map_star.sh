#!/bin/bash

# Required modules
module load fastqc
module load rnastar
module load samtools
module load deeptools

# Initialize options
r1f=""
r2f=""
fname=""
genome=""
outpath=""
threads=1

# Collect options
while getopts "1:2:n:g:c:o:" option; do
	case $option in
		1) r1f="$OPTARG";;
		2) r2f="$OPTARG";;
		n) fname="$OPTARG";;
		g) genome="$OPTARG";;
		c) threads=$OPTARG;;
		o) outpath="$OPTARG";;
		\?) # Invalid option
			echo "Error: Invalid option"
			exit;;
	esac
done

# Check that required options were passed
if [[ -z $r1f ]]; then
	echo "R1 file unspecified"
	exit 22
elif [[ -z $r2f ]]; then
	echo "R2 file unspecified"
	exit 22
elif [[ -z $fname ]]; then
	echo "File name unspecified"
	exit 22
elif [[ -z $genome ]]; then
	echo "Genome unspecified"
	exit 22
elif [[ -z $outpath ]]; then
	echo "Output path unspecified"
	exit 22
fi

# Actual execution
cd $outpath

echo "R1: $r1f"
echo "R2: $r2f"
echo "Name: $fname"
echo "Genome: $genome"
date

mkdir -p 01_fastqc_reports
fastqc -v
fastqc $r1f --outdir=./01_fastqc_reports --threads=$threads
fastqc $r2f --outdir=./01_fastqc_reports --threads=$threads

mkdir -p 02_mapped
echo "STAR " $(STAR --version)
STAR --readFilesIn $r1f $r2f \
	 --genomeDir "$genome" --outFileNamePrefix "./02_mapped/${fname}_" \
	 --runThreadN $threads --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# samtools view -F 0x100 -u - | samtools sort -o "./04_mapped/${fname}.bam" -

# samtools index "./04_mapped/${fname}.bam" > "./04_mapped/${fname}.bam.bai"

# mkdir -p 03_bigwigs
# bamCoverage -b "./04_mapped/${fname}.bam" -o "./06_bigwigs/${fname}.bw" --normalizeUsing CPM -p $threads

date
