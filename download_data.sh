#!/bin/bash
#SBATCH --partition=standard --time=2-00:00:00 -c 6 --mem=10G --output=data_download.log

# This script can take a long time and downloads a lot of data
# It might be preferable to run it in parts

# Peregrine data
wget -nv -P ./data/external_resources/peregrine/ http://data.pantherdb.org/ftp/peregrine_data/PEREGRINEenhancershg38
wget -nv -P ./data/external_resources/peregrine/ http://data.pantherdb.org/ftp/peregrine_data/PEREGRINEtissues
wget -nv -P ./data/external_resources/peregrine/ http://data.pantherdb.org/ftp/peregrine_data/enh_gene_link_assay_tissue_pval_snp_score
wget -nv -P ./data/external_resources/peregrine/ http://data.pantherdb.org/ftp/peregrine_data/assaytable.txt

module load sratoolkit

# Our data, ATAC
prefetch --output-directory ./data/fastq --option-file ./data/fastq/SRR_Acc_List.txt
fasterq-dump --outdir ./data/fastq ./data/fastq/SRR*

# Ucar data, ATAC
# Same as above, but requires access permissions

# Ucar data, RNA
# Same as above, but requires access permissions

# Covid data, ATAC
# Data available on GEO: GSE206266

# Then I compressed with a separate script

# Illumina probe data
wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv -O ./data/external_resources/methylation/full_probe_info_450.csv
cat ./data/external_resources/methylation/full_probe_info_450.csv | tail -n +8 | head -n -851 | cut -d ',' -f 1,11,12,13 > ./data/external_resources/methylation/short_probe_info_450.csv

# Hannum methylation data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz -O ./data/external_resources/methylation/hannum_series_matrix.txt.gz
gzip -d ./data/external_resources/methylation/hannum_series_matrix.txt.gz

# Hannum clock coefficients
https://ars.els-cdn.com/content/image/1-s2.0-S1097276512008933-mmc2.xlsx -> convert first sheet to csv from excel -> hannum_coef.csv

# Horvath clock coefficients
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv -> remove 2 lines -> horvath_coef.csv