The order of analysis is as follows:
- Most of the data used in the paper is collected by ./download_data.sh collects, except for datasets used for validation
- FastQ files for ATAC and RNA seq are aligned and processed using the scripts in pipeline_atac and pipeline_rna
- QC, normalization, corrections are performed by R scripts named clean_data_(atac/re/rna) in ./scripts
- OCR-gene link tables are prepared in ./scripts/peregrine_link_enhancers.R
- Methylation data preprocessing is done in ./scripts/methylation_prepro.R
- Datasets are prepared for clock training inside jupyter notebooks and ncv is run in parallel on slurm
- The core of the analysis and figures are done in R by ./scripts/analysis_analysis.R
- Integrative analysis of ATAC, RNA and MET is done in ./scripts/analysis_integ.R
- Validation of the integrative analysis using a different RNA-seq dataset is done in ./scripts/analysis_heatmap_val.R
- Validation of clock performance on externa data is done in ./scripts/analysis_atac_val.R

Main ATAC pipeline (./pipeline_atac):
- 01_dispatcher_mapping.sh: submits separate slurm jobs for ATAC mapping job, calls map_and_clean.sh
- 02_peakset_and_counts.sh: creates peaksets and count tables
- 03_rep_element_counts.sh: counts reads over repetitive elements
- 04_dispatcher_bigwigs.sh: submits separate slurm jobs for ATAC coverage tracks, calls make_bigwig.sh
- map_and_clean.sh: runs atac pipeline for one sample, called by 01_dispatched_mapping.sh
- make_bigwig.sh: generates pileup bigwigs with appropriate normalization factors

Additional ATAC utilities (./pipeline_atac):
- frag_length_hist.sh: generates histograms of fragment lengths for ATAC reads
- tss_enrichment.sh: creates TSS enrichment profiles
	
Main RNA pipeline (./pipeline_rna):
- 01_dispatcher.sh: submits separate slurm jobs for RNA mapping jobs, calls map_star.sh
- 02_count_reads.sh: creates count table
- map_star: maps RNA reads using STAR for one sample, called by 01_dispatcher.sh

Clocks (./clocks):
- clock_*.ipynb: jupyter notebook used to prepare datasets for ncv
- clock_utils.py: contains various clock scripts. Not all are used
- collect_ncv_results.py: summarized ncv outcome
- final_clock_cv.py: runs a final clock on all data
- ncv_fold.py: runs one fold of ncv
- ncv_fold_with_correction.py: runs one fold of ncv performing cell composition correction within ncv
- run_clock.sh: script used to run ncv in parallel (and train a final clock if desired)

Scripts (./scripts):
- analysis_atac_val.R: validates clock performance on Marquez and Giroux data
- analysis_heatmap_val.R: validates integrative analysis of RNA, ATAC and MET
- analysis_integ.R: main integrative analysis
- analysis_main.R: main analysis
- clean_data_atac.R: perform QC, normalize, correct ATAC data
- clean_data_re.R: perform QC, normalize, collapse ATAC data for repetitive elements
- clean_data_rna.R: perform QC, normalize RNA data
- methylation_prepro.R: calculate age correlation of MET data, match CpGs to OCRs
- peregrine_link_enhancers.R: save gene-OCR links made with PEREGRINE or ChipSeeker
- utils.R: various utility scripts