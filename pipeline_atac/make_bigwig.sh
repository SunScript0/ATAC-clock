#!/bin/bash

# Required modules
module load samtools
module load bedtools
module load deeptools

# Initialize options
f=""
genome=""
outpath=""
qc_summary=""

# Collect options
while getopts "f:g:o:q:" option; do
	case $option in
		f) f="$OPTARG";;
		g) genome="$OPTARG";;
		o) outpath="$OPTARG";;
		q) qc_summary="$OPTARG";;
		\?) # Invalid option
			echo "Error: Invalid option"
			exit;;
	esac
done

# Check that required options were passed
if [[ -z $f ]]; then
	echo "R1 file unspecified"
	exit 22
elif [[ -z $genome ]]; then
	echo "Genome unspecified"
	exit 22
elif [[ -z $outpath ]]; then
	echo "Output path unspecified"
	exit 22
elif [[ -z $qc_summary ]]; then
	echo "QC summary unspecified"
	exit 22
fi

# Actual execution
cd $outpath

echo "File: $f"
echo "Genome: $genome"
date

mkdir -p 08_bigwigs
sname="$(basename "$f")"
samtools sort -n -O BAM "./05_clean_BAMs/${sname}_clean.bam" > "./08_bigwigs/${sname}_qsorted.bam"
bedtools bamtobed -bedpe -i "./08_bigwigs/${sname}_qsorted.bam" > "./08_bigwigs/${sname}.bedpe"
cut -f1,2,6,7 "./08_bigwigs/${sname}.bedpe" |\
	bedtools flank -b 1 -g $genome |\
	bedtools slop -b 100 -g $genome |\
	bedtools bedtobam -g $genome |\
	samtools sort > "./08_bigwigs/${sname}_cut_sites.bam"
samtools index "./08_bigwigs/${sname}_cut_sites.bam"
# # I scale directly to reads in peak x 10e6 / bin_size: (1'000'000 * n) / (rip * bin_size) = (100'000 * n) / rip if bin_size = 10
# # Actualy i multiply by 100 too to get nicer numbers
reads_in_peak=$(grep "$sname" $qc_summary | cut -f9)
bw_scale_factor=$(python -c "print(10000000 / $reads_in_peak)")
echo "Scale factor: $bw_scale_factor"
bamCoverage -bs 10 --scaleFactor $bw_scale_factor -b "./08_bigwigs/${sname}_cut_sites.bam" -o "./08_bigwigs/${sname}.bw"
rm "./08_bigwigs/${sname}_qsorted.bam"
rm "./08_bigwigs/${sname}.bedpe"

date