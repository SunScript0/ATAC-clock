# This script prepares public methylation data for analysis
# - Reads probe info from illumina 450 methylation array
# - Filters the probe info, splits it into coords refering to h36 and h37 genomes
# - Writes the files to use in LiftOver to convert to h38
# - Reads them back and maps them to called peaks#
# - Imports hannum methylation data, computes age correlations and saves them

library(data.table)
library(tidyverse)
library(GenomicRanges)
setwd("/scratch/fmorandi/ChromAcc-clock")
source("./scripts/utils.R")

paths = list()
paths$met = "./data/external_resources/methylation/"
paths$data = "./data/paper_data/"

##### PREP PROBE BED FILES #####

probe_info = read.csv(paste0(paths$met, "short_probe_info_450.csv")) %>%
  drop_na(Genome_Build) %>% # Few sites have no coordinate at all, will ignore them
  mutate(CHR = paste0("chr", CHR))
probes37 = probe_info[probe_info$Genome_Build == 37, c("CHR", "MAPINFO", "MAPINFO", "IlmnID")]
probes36 = probe_info[probe_info$Genome_Build == 36, c("CHR", "MAPINFO", "MAPINFO", "IlmnID")]
# Write tables for liftover
write.table(probes37, file=paste0(paths$met, "probes450_h37.bed"), quote = F, row.names = F, col.names = F)
write.table(probes36, file=paste0(paths$met, "probes450_h36.bed"), quote = F, row.names = F, col.names = F)
# Remove big tables
rm(probe_info, probes36, probes37)

##### LIFTOVER #####

# Go to LiftOver website or use command line tool to convert to h38 genome version

##### MAP CPGS TO OCRS #####

# Then load back into R
probes36_38 = read.table(paste0(paths$met, "probes450_h36_h38.bed"))
probes37_38 = read.table(paste0(paths$met, "probes450_h37_h38.bed"))
probes38 = rbind(probes37_38, probes36_38)[1:4]
colnames(probes38) = c("chr", "start", "end", "probe_id")

# Read atac peaks from processed dataset
peak_info = read.table(paste0(paths$data, "peak_info.tsv"))
peak_ranges = GRanges(
  seqnames = peak_info$Chr,
  ranges = IRanges(peak_info$Start, end = peak_info$End, names = rownames(peak_info))
)
probe_ranges = GRanges(
  seqnames = probes38$chr,
  ranges = IRanges(probes38$start, names = probes38$probe_id)
)

mappings = data.frame(findOverlaps(peak_ranges, probe_ranges))
mappings$queryHits=names(peak_ranges)[mappings$queryHits]
mappings$subjectHits=names(probe_ranges)[mappings$subjectHits]
colnames(mappings) = c("peak_id", "probe_id")
mappings = merge(mappings, probes38[c("chr", "start", "probe_id")], by="probe_id")
colnames(mappings) = c("probe_id", "peak_id", "probe_chr", "probe_coord38")

write.table(mappings, file=paste0(paths$data, "peak_cpg_map.tsv"), row.names = F, quote=F)

rm(mappings, peak_info, peak_ranges, probe_ranges, probes36_38, probes37_38, probes38)

##### PREPROCESS HANNUM #####

# Get first 100 lines of series matrix
header = readLines(paste0(paths$met, "hannum_series_matrix.txt"), n=100)

# Find line with GSM ids
ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

# Find lines with metadata
hannum_meta = header[grepl("!Sample_characteristics", header)]
hannum_meta = gsub('"', "", hannum_meta)
hannum_meta = str_split(hannum_meta, pattern = "\t", simplify = T)
hannum_meta = data.frame(t(hannum_meta))
hannum_meta = hannum_meta[-c(1), ]
colnames(hannum_meta) = c("age", "source", "plate", "sex", "ethnicity", "tissue")
hannum_meta = hannum_meta %>%
  dplyr::select(age, sex) %>%
  mutate(age = str_replace_all(age, "age \\(y\\): ", "")) %>%
  mutate(sex = str_replace_all(sex, "gender: ", "")) %>%
  mutate(age = as.numeric(age))
rownames(hannum_meta) = ids

# Find start of count table
mat_start = which(grepl("series_matrix_table_begin", header))
# Read beta values
hannum_data = fread(paste0(paths$met, "hannum_series_matrix.txt"), skip = mat_start, data.table=F)
rownames(hannum_data) = hannum_data$ID_REF
hannum_data = hannum_data[, -c(1)]
hannum_data = data.frame(t(hannum_data))

# Calculate spearman correlations and pvalues
hannum_age_cors = cor_test_matrix(hannum_data, hannum_meta$age, "spearman")
hannum_age_cors$padj = p.adjust(hannum_age_cors$pvals, method = "fdr")
colnames(hannum_age_cors) = c("corr","pval","padj")

write.table(hannum_age_cors, paste0(paths$data, "hannum_age_cors.tsv"), quote=F)
