#!/bin/bash

# Required modules
module load fastqc
module load trimgalore
module load cutadapt
module load bowtie2
module load picard
module load samtools
module load bedtools

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

# mkdir -p 01_raw_fastqc_reports
# fastqc -v
# fastqc $r1f --outdir=./01_raw_fastqc_reports --threads=$threads
# fastqc $r2f --outdir=./01_raw_fastqc_reports --threads=$threads

# mkdir -p 02_trimmed_fastq
# echo "TrimGalore version:"
# trim_galore -v
# trim_galore --output_dir ./02_trimmed_fastq --paired $r1f $r2f --basename "${fname}"

# # Stop pipeline if fastq is corrupted
# if test -f "./02_trimmed_fastq/${fname}_R1_val_1.fq.gz"; then
    # echo "Fastqs trimmed succesfully: continuing"
# else 
    # echo "Fastqs failed trimming, likely corrupted: exiting"
	# date
	# exit
# fi

# mkdir -p 03_trimmed_fastqc_reports
# fastqc "./02_trimmed_fastq/${fname}_R1_val_1.fq.gz" --outdir=./03_trimmed_fastqc_reports --threads=$threads
# fastqc "./02_trimmed_fastq/${fname}_R2_val_2.fq.gz" --outdir=./03_trimmed_fastqc_reports --threads=$threads

# mkdir -p 04_mapped
# echo "Bowtie2 " $(bowtie2 --version | head -n 1)
# bowtie2 --very-sensitive -X 1000 --dovetail -p $threads \
# -x $genome \
# -1 "./02_trimmed_fastq/${fname}_R1_val_1.fq.gz" \
# -2 "./02_trimmed_fastq/${fname}_R2_val_2.fq.gz" \
# | samtools view -u - | samtools sort -o "./04_mapped/${fname}.bam" -

mkdir -p 05_clean_BAMs
# Print flagstat to get number of raw pairs, mapped pairs, proper pairs
samtools flagstat "./04_mapped/${fname}.bam"
# Remove multimapping and keep only proper pairs
samtools view -f 0x2 -F 0x100 -b "./04_mapped/${fname}.bam" > "./05_clean_BAMs/${fname}_propairs.bam"
# Make bam index
samtools index "./05_clean_BAMs/${fname}_propairs.bam" > "./05_clean_BAMs/${fname}_propairs.bam.bai"
# Print number of mitochondrial pairs from bam index
echo "Found $(($(samtools idxstats "./05_clean_BAMs/${fname}_propairs.bam" | grep "chrM" | cut -f 3) / 2)) mitochondrial pairs"
# Duplicate metrics and remove all duplicates
java -jar $PICARD MarkDuplicates I="./05_clean_BAMs/${fname}_propairs.bam" \
	O="./05_clean_BAMs/${fname}_dedup.bam" \
	M="./05_clean_BAMs/${fname}_duplicate_metrics.txt" \
	REMOVE_DUPLICATES=true
# Remove propairs bam and index
rm "./05_clean_BAMs/${fname}_propairs.bam"
rm "./05_clean_BAMs/${fname}_propairs.bam.bai"
# Index dedup and remove chrM reads
samtools index "./05_clean_BAMs/${fname}_dedup.bam" > "./05_clean_BAMs/${fname}_dedup.bam.bai"
samtools idxstats "./05_clean_BAMs/${fname}_dedup.bam" | cut -f 1 | grep -v chrM | xargs samtools view -b "./05_clean_BAMs/${fname}_dedup.bam" > "./05_clean_BAMs/${fname}_clean.bam"
samtools index "./05_clean_BAMs/${fname}_clean.bam" > "./05_clean_BAMs/${fname}_clean.bam.bai"
# Remove no dupes bam and index
rm "./05_clean_BAMs/${fname}_dedup.bam"
rm "./05_clean_BAMs/${fname}_dedup.bam.bai"
# Print number of clear pairs
echo "Clean BAM contains $(($(samtools view -c "./05_clean_BAMs/${fname}_clean.bam") / 2)) pairs"

mkdir -p 06_macs2_outputs
bedtools bamtobed -i "./05_clean_BAMs/${fname}_clean.bam" | gzip > "./06_macs2_outputs/${fname}.bed.gz"

module load macs
macs2 callpeak -t "./06_macs2_outputs/${fname}.bed.gz" \
	--outdir "./06_macs2_outputs" --name "$fname" \
	-f BED -g "hs" --keep-dup "all" -q 0.01 \
	--nomodel --shift -100 --extsize 200

date