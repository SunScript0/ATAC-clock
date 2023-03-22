library(data.table)
library(tidyverse)

setwd("/scratch/fmorandi/ChromAcc-clock")

paths = list()
paths$pipeline ="./data/pipeline_outputs/atac/"
paths$out="./data/paper_data/"

##### LOAD DATA #####

# RE only table
counts_res = fread(paste0(paths$pipeline, "07_peaksets_and_tables/counts_rtes.tsv"))
counts_res = counts_res %>%
  mutate(fam = str_extract(counts_res$Geneid, "^(\\w*[^/]?)\\/", group=1)) %>%
  group_by(fam) %>%
  dplyr::select(-c(Geneid:Length)) %>%
  summarise_all(sum)
counts_res = data.frame(counts_res)
rownames(counts_res) = counts_res$fam
counts_res = counts_res[, -1]
colnames(counts_res) = str_extract(colnames(counts_res), "(SRR[0-9]+)", group=1)
counts_res = data.frame(t(counts_res))

# RE only featureCounts report
fc_report_re = read.table(paste0(paths$pipeline, "07_peaksets_and_tables/counts_rtes.tsv.summary"),
                       header=T, row.names=1)
colnames(fc_report_re) = str_extract_all(colnames(fc_report_re), "(SRR[0-9]+)")
fc_report_re = data.frame(t(fc_report_re))
fc_report_re$Total = rowSums(fc_report_re)

# RE vs OCR table
counts_res_ocr = fread(paste0(paths$pipeline, "07_peaksets_and_tables/counts_ocrs_and_rtes.tsv"))
counts_res_ocr = counts_res_ocr %>%
  mutate(RE = grepl("/", Geneid)) %>%
  dplyr::select(-c(Geneid:Length)) %>%
  group_by(RE) %>%
  summarise_all(sum)
counts_res_ocr = data.frame(counts_res_ocr)
rownames(counts_res_ocr) = c("OCR", "Repetitive")
counts_res_ocr = counts_res_ocr[, -1]
colnames(counts_res_ocr) = str_extract(colnames(counts_res_ocr), "(SRR[0-9]+)", group=1)
counts_res_ocr = data.frame(t(counts_res_ocr))

# RE vs OCR featureCounts report
fc_report_re_ocr = read.table(paste0(paths$pipeline, "07_peaksets_and_tables/counts_ocrs_and_rtes.tsv.summary"),
                          header=T, row.names=1)
colnames(fc_report_re_ocr) = str_extract_all(colnames(fc_report_re_ocr), "(SRR[0-9]+)")
fc_report_re_ocr = data.frame(t(fc_report_re_ocr))
fc_report_re_ocr$Total = rowSums(fc_report_re_ocr)

meta = read.table(paste0(paths$out, "meta_final.tsv"), header=T)

##### RENAME SAMPLES #####

dict = data.frame(Subject = meta$Subject, SRR = meta$SRR_atac)
dict = drop_na(dict)
rownames(dict) = dict$SRR

rownames(counts_res) = dict[rownames(counts_res), "Subject"]
rownames(counts_res_ocr) = dict[rownames(counts_res_ocr), "Subject"]
rownames(fc_report_re) = dict[rownames(fc_report_re), "Subject"]
rownames(fc_report_re_ocr) = dict[rownames(fc_report_re_ocr), "Subject"]

##### NORMALIZE #####

cpm_re = 1e6 * sweep(counts_res, 1, fc_report_re$Total, "/")
cpm_re_ocr = 100 * sweep(counts_res_ocr, 1, fc_report_re_ocr$Total, "/")

##### SAVE OUTPUTS #####

write.table(cpm_re, file=paste0(paths$out, "data_atac_re_cpm.tsv"), quote=F, sep="\t")
write.table(cpm_re_ocr, file=paste0(paths$out, "data_atac_re_ocr_perc.tsv"), quote=F, sep="\t")
