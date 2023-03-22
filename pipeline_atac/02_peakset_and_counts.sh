#!/bin/bash
#SBATCH --partition=standard --time=24:00:00 -c 4 --mem=100G --output=covid_counts.log

module load bedtools
module load subread

blacklist=/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/hg38-blacklist.v2.bed
basepath=/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/GSE206266_atac
precomputed_peakset=/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac/07_peaksets_and_tables/peaks_filt.saf
threads=4

cd $basepath

date

mkdir -p 07_peaksets_and_tables/
# # Concatenate all narrowPeak files and sort
# cat ./06_macs2_outputs/*_peaks.narrowPeak | cut -f1,2,3,4 | bedtools sort > ./07_peaksets_and_tables/peaks_cat.bed
# # Bedtools merge on concatenated peak file
# bedtools merge -i ./07_peaksets_and_tables/peaks_cat.bed > ./07_peaksets_and_tables/peaks_merged.bed

# # Multiinter on all narrowPeak files
# bedtools multiinter -i ./06_macs2_outputs/*.narrowPeak > ./07_peaksets_and_tables/multinter.txt
# # Select regions which were called as peaks in at least  samples
# cat ./07_peaksets_and_tables/multinter.txt | awk '$4>50' > ./07_peaksets_and_tables/multinter_thresholded.txt

# # Only keep peaks which overlap regions called as peak in at least 50 samples
# bedtools intersect -wa -u -sorted \
	# -a ./07_peaksets_and_tables/peaks_merged.bed \
	# -b ./07_peaksets_and_tables/multinter_thresholded.txt \
	# > ./07_peaksets_and_tables/peaks_merged_filt.bed

# # Remove blacklist
# bedtools intersect -v \
	# -a ./07_peaksets_and_tables/peaks_merged_filt.bed \
	# -b $blacklist > ./07_peaksets_and_tables/peaks_filt.bed
	
# # Make SAF for featureCounts
# awk 'OFS="\t" {print $1":"$2"-"$3, $1, $2, $3, "."}' ./07_peaksets_and_tables/peaks_filt.bed > ./07_peaksets_and_tables/peaks_filt.saf
# # Count 5' ends of reads over peaks
# featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	# -a ./07_peaksets_and_tables/peaks_filt.saf \
	# -o ./07_peaksets_and_tables/counts.tsv \
	# ./05_clean_BAMs/*.bam
featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	-a $precomputed_peakset \
	-o ./07_peaksets_and_tables/counts.tsv \
	./05_clean_BAMs/*.bam
# Not counting fractional because 5' end is a single bp and there should be no multiple overlapping peaks
# featureCounts -a $saf -o $out -F "SAF" -O -M --fraction -p --countReadPairs -B $ins

# rm ./07_peaksets_and_tables/peaks_cat.bed
# rm ./07_peaksets_and_tables/multinter.txt
# rm ./07_peaksets_and_tables/multinter_thresholded.txt
# rm ./07_peaksets_and_tables/peaks_merged_filt.bed

date