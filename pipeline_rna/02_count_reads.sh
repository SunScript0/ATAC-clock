#!/bin/bash
#SBATCH --partition=standard --time=6:00:00 -c 2 --mem=80G --output=rna_counts.log

module load subread

basepath=/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/rna_ucar
gtf=/scratch/fmorandi/external/references/GRCh38-hg38-UCSC/hg38.ensGene.gtf
threads=2

cd $basepath

date

mkdir -p 03_counts

featureCounts -F "GTF" -p --countReadPairs -B -T $threads \
	-a $gtf \
	-o ./03_counts/counts.tsv \
	./02_mapped/*.bam
	
date