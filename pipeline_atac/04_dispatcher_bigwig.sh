#!/bin/bash
# Job settings
partition="standard"
max_time="4:00:00"
cores=1
mem="5G"

files="$(find "/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac/05_clean_BAMs" -name "*_clean.bam" | sed 's/_clean.bam//' | sort)"
don="$(find /scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac/08_bigwigs -name *.bw -printf "%f\n")"
don="$(echo "$don"  | sed "s/.bw//" | tr "\n" "|" | sed "s/|/\\\|/g" | sed "s/\\\|$//")"
files="$(echo "$files" | grep -v "$don")"

echo "$files"

# Script settings
script="/scratch/fmorandi/ChromAcc-clock/pipeline_atac/make_bigwig.sh"
basepath=/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac
genome_index=/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/hg38.fa.fai
qc_summary=/scratch/fmorandi/ChromAcc-clock/data/paper_data/qc_summary.tsv

cd $basepath
mkdir -p B_bigwig_logs

for f in $files
do
	fname=$(basename $f)
	logf="${basepath}/B_bigwig_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -f "$f" -g "$genome_index" -q "$qc_summary" -o "$basepath"
done