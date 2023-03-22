#!/bin/bash
#SBATCH --partition=standard --time=6:00:00 -c 2 --mem=80G --output=re_counts.log

module load bedtools
module load subread

basepath=/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac
repmasked=/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/RepeatMaskerOut/hg38.fa.out
threads=2

cd $basepath

date

# # Make SAF from RepeatMaskerOut output
# # Family (superfamily/family/repeat_id, chr, start, end, strand
# tail -n +4 $repmasked | awk 'OFS="\t" {print $11"/"$10"/"$15, $5, $6, $7, "."}' > ./07_peaksets_and_tables/rtes.saf

# # Concatenate peakset and repeat regions
# cat ./07_peaksets_and_tables/peaks_filtered.saf ./07_peaksets_and_tables/rtes.saf > ./07_peaksets_and_tables/peaks_and_rtes.saf

# # Count 5' ends of reads over REs
# featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	# -a ./07_peaksets_and_tables/rtes.saf \
	# -o ./07_peaksets_and_tables/counts_rtes.tsv \
	# ./05_clean_BAMs/*.bam
	
# Count 5' ends of reads over concatenated peaks and REs
# This time I will need to account for overlapping features
featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	-O --fraction \
	-a ./07_peaksets_and_tables/peaks_and_rtes.saf \
	-o ./07_peaksets_and_tables/counts_ocrs_and_rtes.tsv \
	./05_clean_BAMs/*.bam

date