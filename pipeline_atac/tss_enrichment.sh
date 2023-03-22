tss_regions="/scratch/fmorandi/external/references/hg38/hg38_TSS.bed"
bw_path="/scratch/fmorandi/ChromAcc-clock/data/pipeline_outputs/atac/08_bigwigs"

echo "Starting to compute matrix"

computeMatrix reference-point -S "$bw_path/"*.bw \
							  -R "$tss_regions" \
						      -b 2000 -a 2000 -o "g1-2_matrix_all.gz"

echo "Done"
							  
# plotHeatmap -m g1-2_matrix_all.gz -o g1-2_heatmap_all.png
# plotProfile -m g1-2_matrix_all.gz -o g1-2_profile_all.png --outFileNameData g1-2_profile_data_all 