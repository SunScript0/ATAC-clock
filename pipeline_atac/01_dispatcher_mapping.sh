#!/bin/bash
# Job settings
partition="standard"
max_time="32:00:00" # 24h
cores=6 #2
mem="40G" #30G

# File list
files="$(find /scratch/fmorandi/ChromAcc-clock/data/fastq_GSE206266 -name *.fastq.gz)"
files="$(echo "$files" | grep -E "SRR19681105|SRR19681121|SRR19681125|SRR19681134|SRR19681139")"
# files="$(find /scratch/fmorandi/ChromAcc-clock/data/fastq -name SRR17466642*.fastq.gz)"
# don="$(find /scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac/06_macs2_outputs -name *_peaks.narrowPeak -printf "%f\n")"
# don="$(echo "$don"  | sed "s/_peaks.narrowPeak//" | tr "\n" "|" | sed "s/|/\\\|/g" | sed "s/\\\|$//")"
# files="$(echo "$files" | grep -v "$don")"

# echo "$files"
# exit

# Script settings
script="/scratch/fmorandi/ChromAcc-clock/pipeline_atac/map_and_clean.sh"
genome="/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/Bowtie2/GRCh38"
outpath="/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/GSE206266_atac"
r1_pattern="_1\.fastq.gz"
r2_pattern="_2\.fastq.gz"

cd $outpath
mkdir -p A_mapping_logs

r1_files=$(echo "$files" | grep "$r1_pattern")

for r1f in $r1_files
do
	r2f=$(echo $r1f | sed "s/$r1_pattern/$r2_pattern/")
	if [ $(echo $files | grep -c $r2f) -ne 1 ]; then
		echo "R2 file of $r1f not found"
		continue
	fi
	fname=$(basename $r1f | sed "s/$r1_pattern//")
	logf="${outpath}/A_mapping_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -1 "$r1f" -2 "$r2f" -n "$fname" -g "$genome" -c $cores -o "$outpath"
	sleep 0.2
done