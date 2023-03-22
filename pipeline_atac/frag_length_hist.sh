bam_path="/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/05_clean_BAMs"
out_path="/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/extra_qc"

bams="$(find "$bam_path" -name '*_clean.bam' | sort)"
cd $out_path

for bam in $bams
do
	fname="$(basename $bam | sed -r 's/_clean.bam//')"
	bamPEFragmentSize \
	-b "$bam" \
	-o "$fname"_frag_size_hist.png \
	--maxFragmentLength 1000 \
	--samplesLabel "$fname"
done