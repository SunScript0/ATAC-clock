#!/bin/bash
# Job settings
partition="standard"
max_time="5:00:00"
cores=2
mem="40G"

# File list
files="$(find /gpfs/fs2/scratch/fmorandi/ChromAcc-clock/data/fastq_ucar_rna -name *.gz)"
# files="$(echo "$files" | sort | head -n 2)"
# echo "$files"
# exit

# Script settings
script="/scratch/fmorandi/ChromAcc-clock/pipeline_rna/map_star.sh"
genome="/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/STAR_with_SJs"
outpath="/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/rna_ucar"
r1_pattern="_1\.fastq\.gz"
r2_pattern="_2\.fastq\.gz"

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
	logf="./A_mapping_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -1 "$r1f" -2 "$r2f" -n "$fname" -g "$genome" -c $cores -o "$outpath"
done