library(tidyverse)
library(patchwork)
library(scales)
library(sp)

setwd("/scratch/fmorandi/ChromAcc-clock")

source("./scripts/utils.R")

paths = list()
paths$pipeline1 ="./data/pipeline_outputs/atac_ucar/"
paths$pipeline2 ="./data/pipeline_outputs/atac_covid/"
paths$meta1 = "./data/meta_ucar.tsv"
paths$meta2 = "./data/meta_covid2.txt"
paths$out = "./data/paper_addon_atac/"

# Options
th_counts = 1.1e7
th_frip = 0.18
w = 174 # mm
h = 230

dir.create(paste0(paths$out, "figures"), showWarnings = F)

##### PREPRO UCAR DATA #####
###### LOAD DATA ######

counts = read.table(paste0(paths$pipeline1, "07_peaksets_and_tables/counts_our_peaks.tsv"), header=T)
peak_info = counts[, 1:6]
counts = counts[, -c(1:6)]
colnames(counts) = str_extract(colnames(counts), "(SRR[0-9]+)", group=1)
rownames(counts) = paste0("peak", 1:nrow(counts))
rownames(peak_info) = rownames(counts)
counts = data.frame(t(counts))

fc_report = read.table(paste0(paths$pipeline1, "07_peaksets_and_tables/counts_our_peaks.tsv.summary"),
                       header=T, row.names=1)
colnames(fc_report) = str_extract_all(colnames(fc_report), "(SRR[0-9]+)")

meta = read.table(paste0(paths$meta1), header=T, sep="\t", row.names = 1)
meta$Sex = fct_recode(as.factor(meta$Sex), "M" = "1", "F" = "2")
meta = subset(meta, ANALYTE_TYPE == "ATAC-seq")

###### EXTRACT QC INFO #####

logs = dir(paste0(paths$pipeline1, "A_mapping_logs"), full.names=T)

qc = data.frame()
for (log in logs) {
  f = read_file(log)
  sname = str_extract(log, "(SRR[0-9]+)", group=1)
  # TrimGalore! prints the total number of pairs processed after validation
  #   This is the number of raw pairs
  #   Sequences that were trimmed too much are removed
  s = str_extract(f, "Total number of sequences analysed: ([0-9]+)", group=1)
  qc[sname, "raw"] = as.numeric(s)
  # Bowtie2 summary includes the number of pairs processed
  #   This corresponds to the number of pairs after trimming
  s = str_extract(f, "([0-9]+) reads; of these:", group=1)
  qc[sname, "trimmed"] = as.numeric(s)
  # Flagstat prints the number of properly paired reads
  #   I can confirm that this is the number of reads seen from now on
  #   By comparing to the remove_duplicates log
  s = str_extract(f, "([0-9]+) \\+ [0-9]+ properly paired", group=1)
  qc[sname, "proper_pair"] = as.numeric(s) / 2
  # I had the pipeline print the number of mito reads from idxstats
  s = str_extract(f, "Found ([0-9]+) mitochondrial pairs", group=1)
  qc[sname, "mitochondrial"] = as.numeric(s)
  # Picard MarkDuplicates prints the number of dupes in the log (divide by 2 to get pairs)
  s = str_extract(f, "Marking ([0-9]+) records as duplicates", group=1)
  qc[sname, "duplicates"] = as.numeric(s) / 2
  # I had the pipeline print the number of clean reads 
  s = str_extract(f, "Clean BAM contains ([0-9]+) pairs", group=1)
  qc[sname, "clean"] = as.numeric(s)
}

# Verify that qc and fc_report match up
fc_report = fc_report[, rownames(qc)]
all(colnames(fc_report) == rownames(qc))
all(colSums(fc_report)/2 == qc$clean)
# Get frip
fc_report = data.frame(t(fc_report))
qc$cut_sites = fc_report$Assigned + fc_report$Unassigned_NoFeatures
qc$cut_sites_in_peak = fc_report$Assigned
qc$frip = qc$cut_sites_in_peak / qc$cut_sites

###### REMOVE LOW DEPTH SAMPLES #####

ggplot(qc, aes(x=clean, y=frip))+
  geom_point()+
  geom_vline(xintercept=th_counts, color="red")+
  geom_hline(yintercept=th_frip, color="red")+
  labs(x="N of filtered pairs", y="FRIP") +
  theme_light()
ggsave(paste0(paths$out, "figures/qc_atac_thresholds_ucar.pdf"), height = 0.3*h, width = 0.5*w, units="mm", dpi = 300)

# For our data we had filtered by both frip and counts but here, reads are 
# counted on our peaks so FRIP will have been estimated differently.
# So we just filter our samples with < 1.1e7 clean reads

passing = rownames(qc)[qc$clean>1.1e7]
meta = subset(meta, SRR %in% passing)
counts = counts[passing, ]

###### COLLAPSE REPLICATES #####

rownames(meta) = meta$SRR
meta = meta[rownames(counts), ]
all(rownames(counts) == rownames(meta))

groups = meta$Subject
collapsed = data.frame()
for (g in unique(groups)) {
  inds = which(groups == g)
  c_row = colSums(counts[inds, ])
  collapsed[as.character(g), names(c_row)] = c_row
}

rownames(collapsed) = paste0("S", rownames(collapsed))
meta = meta %>%
  mutate(Subject = paste0("S", Subject)) %>%
  dplyr::select(Subject, Sex, Age) %>%
  distinct()
rownames(meta) = NULL

###### TPM NORMALIZATION #####

tpm = 1000 * sweep(collapsed, 2, peak_info$Length, "/")
tpm = 1e6 * sweep(tpm, 1, rowSums(tpm), "/")

###### MARK OUTLIERS #####

pca = prcomp(log1p(tpm))
pca = merge(pca$x, meta, by.x=0, by.y="Subject")
p = ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point()+
  stat_ellipse(level=0.99)
ellipse = ggplot_build(p)$data[[2]]
pca$outlier = as.logical(!point.in.polygon(pca$PC1, pca$PC2, ellipse$x, ellipse$y))
ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point(aes(color=outlier))+
  stat_ellipse(level=0.99)+
  theme_light()
ggsave(paste0(paths$out, "figures/qc_atac_ellipse_ucar.pdf"), height = 0.3*h, width = 0.6*w, units="mm", dpi = 300)

meta$Outlier = pca$outlier

###### SAVE OUTPUTS #####

write.table(tpm, file=paste0(paths$out, "atac_tpm_ucar.tsv"), quote=F, sep="\t")

write.table(peak_info, file=paste0(paths$out, "peak_info.tsv"), quote=F, sep="\t")
write.table(meta, file=paste0(paths$out, "meta_ucar.tsv"), quote=F, sep="\t")
write.table(qc, file=paste0(paths$out, "qc_summary_ucar.tsv"), quote=F, sep="\t")

rm(collapsed, counts, ellipse, fc_report, meta, p, pca, qc, tpm)

##### PREPRO COVID DATA #####
###### LOAD DATA ######

peak_info = read.table(paste0(paths$out, "peak_info.tsv"))
counts = read.table(paste0(paths$pipeline2, "07_peaksets_and_tables/counts.tsv"), header=T)
counts = counts[, -c(1:6)]
colnames(counts) = str_extract(colnames(counts), "(SRR[0-9]+)", group=1)
rownames(counts) = rownames(peak_info)
counts = data.frame(t(counts))

fc_report = read.table(paste0(paths$pipeline2, "07_peaksets_and_tables/counts.tsv.summary"),
                       header=T, row.names=1)
colnames(fc_report) = str_extract_all(colnames(fc_report), "(SRR[0-9]+)")

meta = read.table(paste0(paths$meta2), header=T, sep="\t")
meta = meta[1:10]

###### EXTRACT QC INFO #####

logs = dir(paste0(paths$pipeline2, "A_mapping_logs"), full.names=T)

qc = data.frame()
for (log in logs) {
  f = read_file(log)
  sname = str_extract(log, "(SRR[0-9]+)", group=1)
  # TrimGalore! prints the total number of pairs processed after validation
  #   This is the number of raw pairs
  #   Sequences that were trimmed too much are removed
  s = str_extract(f, "Total number of sequences analysed: ([0-9]+)", group=1)
  qc[sname, "raw"] = as.numeric(s)
  # Bowtie2 summary includes the number of pairs processed
  #   This corresponds to the number of pairs after trimming
  s = str_extract(f, "([0-9]+) reads; of these:", group=1)
  qc[sname, "trimmed"] = as.numeric(s)
  # Flagstat prints the number of properly paired reads
  #   I can confirm that this is the number of reads seen from now on
  #   By comparing to the remove_duplicates log
  s = str_extract(f, "([0-9]+) \\+ [0-9]+ properly paired", group=1)
  qc[sname, "proper_pair"] = as.numeric(s) / 2
  # I had the pipeline print the number of mito reads from idxstats
  s = str_extract(f, "Found ([0-9]+) mitochondrial pairs", group=1)
  qc[sname, "mitochondrial"] = as.numeric(s)
  # Picard MarkDuplicates prints the number of dupes in the log (divide by 2 to get pairs)
  s = str_extract(f, "Marking ([0-9]+) records as duplicates", group=1)
  qc[sname, "duplicates"] = as.numeric(s) / 2
  # I had the pipeline print the number of clean reads 
  s = str_extract(f, "Clean BAM contains ([0-9]+) pairs", group=1)
  qc[sname, "clean"] = as.numeric(s)
}

# Verify that qc and fc_report match up
fc_report = fc_report[, rownames(qc)]
all(colnames(fc_report) == rownames(qc))
all(colSums(fc_report)/2 == qc$clean)
# Get frip
fc_report = data.frame(t(fc_report))
qc$cut_sites = fc_report$Assigned + fc_report$Unassigned_NoFeatures
qc$cut_sites_in_peak = fc_report$Assigned
qc$frip = qc$cut_sites_in_peak / qc$cut_sites

###### REMOVE LOW DEPTH SAMPLES #####

ggplot(qc, aes(x=clean, y=frip))+
  geom_point()+
  geom_vline(xintercept=th_counts, color="red")+
  geom_hline(yintercept=th_frip, color="red")+
  labs(x="N of filtered pairs", y="FRIP") +
  theme_light()
ggsave(paste0(paths$out, "figures/qc_atac_thresholds_covid.pdf"), height = 0.3*h, width = 0.5*w, units="mm", dpi = 300)

# Here again, no filtering by frip and no samples are below the depth threshold

###### COLLAPSE REPLICATES #####

rownames(meta) = meta$Run
meta = meta[rownames(counts), ]
all(rownames(counts) == rownames(meta))

groups = paste(meta$Subject, meta$Timepoint, sep="_T")
collapsed = data.frame()
for (g in unique(groups)) {
  inds = which(groups == g)
  c_row = colSums(counts[inds, ])
  collapsed[as.character(g), names(c_row)] = c_row
}

meta = meta %>%
  dplyr::select(-Rep, -Run) %>%
  distinct() %>%
  mutate(Id = paste(Subject, Timepoint, sep="_T"), .before=Subject)
rownames(meta) = NULL

###### TPM NORMALIZATION #####

tpm = 1000 * sweep(collapsed, 2, peak_info$Length, "/")
tpm = 1e6 * sweep(tpm, 1, rowSums(tpm), "/")

###### MARK OUTLIERS #####

tpm_log = log1p(tpm)
pca = prcomp(tpm_log)
pca = merge(pca$x, meta, by.x=0, by.y="Id")
p = ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point()+
  stat_ellipse(level=0.99)
ellipse = ggplot_build(p)$data[[2]]
pca$outlier = as.logical(!point.in.polygon(pca$PC1, pca$PC2, ellipse$x, ellipse$y))
ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point(aes(color=outlier))+
  stat_ellipse(level=0.99)+
  theme_light()
ggsave(paste0(paths$out, "figures/qc_atac_ellipse_covid.pdf"), height = 0.3*h, width = 0.6*w, units="mm", dpi = 300)

# meta$Outlier = pca$outlier
# No outliers

###### SAVE OUTPUTS #####

write.table(tpm, file=paste0(paths$out, "atac_tpm_covid.tsv"), quote=F, sep="\t")

write.table(meta, file=paste0(paths$out, "meta_covid.tsv"), quote=F, sep="\t")
write.table(qc, file=paste0(paths$out, "qc_summary_covid.tsv"), quote=F, sep="\t")

##### CHECKPOINT #####

library(tidyverse)
library(data.table)
library(ggpubr)
library(ggridges)
library(patchwork)
setwd("/scratch/fmorandi/ChromAcc-clock")
source("./scripts/utils.R")

paths = list()
paths$data_val = "./data/paper_addon_atac/"
paths$data_us = "./data/paper_data/"
paths$clocks_main1 = "./clocks/parallel/2023-02-03_14-38_tpm/"
paths$clocks_main2 = "./clocks/parallel/2023-03-19_10-06_tpm_all_samples/"
paths$out = "./data/paper_data/outputs_2023-03-20_13-43-53/"

# Our ATAC data
atac_us = fread(paste0(paths$data_us, "data_atac_tpm.tsv"))
atac_us = data.frame(atac_us, row.names = 1)
atac_us_log = log1p(atac_us)
# Ucar ATAC data
atac_ucar = fread(paste0(paths$data_val, "atac_tpm_ucar.tsv"))
atac_ucar = data.frame(atac_ucar, row.names = 1)
atac_ucar_log = log1p(atac_ucar)
# Covid ATAC data
atac_covid = fread(paste0(paths$data_val, "atac_tpm_covid.tsv"))
atac_covid = data.frame(atac_covid, row.names = 1)
atac_covid_log = log1p(atac_covid)
rm(atac_us, atac_ucar, atac_covid)

# Our meta
meta_us = read.table(paste0(paths$data_us, "meta_final.tsv"))
rownames(meta_us) = meta_us$Subject
# Ucar meta
meta_ucar = read.table(paste0(paths$data_val, "meta_ucar.tsv"))
rownames(meta_ucar) = meta_ucar$Subject
# Covid meta
meta_covid = read.table(paste0(paths$data_val, "meta_covid.tsv"))
rownames(meta_covid) = meta_covid$Id

# Load clock 1 coefficients
final_coefs1 = fread(paste0(paths$clocks_main1, "final_coefs.tsv"), header = T)
final_coefs1 = data.frame(final_coefs1, row.names = 1)
clock_intercept1 = final_coefs1["intercept", "coef"]
final_coefs1 = final_coefs1[final_coefs1$coef != 0, ]
final_coefs1 = final_coefs1[-which(rownames(final_coefs1) == "intercept"), ]
clock_peaks1 = rownames(final_coefs1)
# Load clock 2 coefficients
final_coefs2 = fread(paste0(paths$clocks_main2, "final_coefs.tsv"), header = T)
final_coefs2 = data.frame(final_coefs2, row.names = 1)
clock_intercept2 = final_coefs2["intercept", "coef"]
final_coefs2 = final_coefs2[final_coefs2$coef != 0, ]
final_coefs2 = final_coefs2[-which(rownames(final_coefs2) == "intercept"), ]
clock_peaks2 = rownames(final_coefs2)

# Plotting settings
theme_set(theme_light(base_size=9))
theme_update(panel.grid.major = element_blank())
theme_update(panel.grid.minor = element_blank())
theme_update(plot.title = element_text(hjust = 0.5, size=9))
theme_update(axis.text = element_text(size = 9, color="black"))
update_geom_defaults("text", list(size = 7*0.35))
update_geom_defaults("point", list(size = 1))
# A4 is 210 x 297
w = 210 - 25.4*2 # was 174
h = 297 - 25.4*2 # was 230
# Colors
cols = list()
cols$azure = "#0171BD"
cols$light_blue = "#95C1E3"
cols$blue = "#0248A2"
cols$light_red = "#D25C5F"
cols$red = "#8A0605"
cols$dark_red = "#550D0A"

##### TEST FINAL CLOCK #####

###### TEST ON US ######
tmp = atac_us_log[, clock_peaks1]
tmp = sweep(tmp, 2, final_coefs1$mu, "-")
tmp = sweep(tmp, 2, sqrt(final_coefs1$var), "/")
preds = rowSums(sweep(tmp, 2, final_coefs1$coef, "*")) + clock_intercept1
meta_us$Preds = preds[meta_us$Subject]
ggplot(meta_us, aes(x=Age, y=Preds, color=PassesQC_atac))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  guides(color="none")
ggsave(paste0(paths$out, "clock_val_us.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

###### TEST ON UCAR ###### 
tmp = atac_ucar_log[, clock_peaks1]
tmp = sweep(tmp, 2, final_coefs1$mu, "-")
tmp = sweep(tmp, 2, sqrt(final_coefs1$var), "/")
preds = rowSums(sweep(tmp, 2, final_coefs1$coef, "*")) + clock_intercept1
meta_ucar$Preds = preds[meta_ucar$Subject]
meta_ucar$Error = meta_ucar$Preds - meta_ucar$Age
anno_text = sprintf("RMSE = %.2f\nMAE = %.2f\nr = %.2f", 
                    sqrt(mean(meta_ucar$Error^2)),
                    median(abs(meta_ucar$Error)),
                    cor(meta_ucar$Age, meta_ucar$Preds))
lims = c(min(c(meta_ucar$Age, meta_ucar$Preds)),
         max(c(meta_ucar$Age, meta_ucar$Preds)))
ggplot(meta_ucar, aes(x=Age, y=Preds))+
  geom_point(color=cols$azure)+
  geom_abline(slope=1, color=cols$dark_red)+
  labs(y="Predicted Age")+
  annotate("text", x = 23, y=100, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  lims(x=lims, y=lims)
ggsave(paste0(paths$out, "clock_val_ucar.pdf"), height = 60.7, width = w/2, units="mm", dpi = 300)

###### TEST ON COVID ###### 
tmp = atac_covid_log[, clock_peaks2]
tmp = sweep(tmp, 2, final_coefs2$mu, "-")
tmp = sweep(tmp, 2, sqrt(final_coefs2$var), "/")
preds = rowSums(sweep(tmp, 2, final_coefs2$coef, "*")) + clock_intercept2
meta_covid$Preds = preds[meta_covid$Id]
meta_covid$Error = meta_covid$Preds - meta_covid$Age
model = lm(Error~Age, data=meta_covid)
meta_covid$Error_corrected = meta_covid$Error - meta_covid$Age * model$coefficients[2] - model$coefficients[1]

ggplot(meta_covid, aes(x=Age, y=Preds))+
  geom_point(color=cols$azure)+
  geom_abline(slope=1, color=cols$dark_red)+
  labs(y="Predicted Age")+
  stat_cor()+
  lims(x=c(10, 80), y=c(10, 80))
ggsave(paste0(paths$out, "clock_val_covid.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

meta_covid2 = meta_covid %>%
  mutate(PCR = fct_recode(PCR, "Negative" = "Negative*")) %>%
  mutate(Affected = Group %in% c("Mild", "Moderate"))

ggplot(meta_covid2, aes(x=PCR, y=Error, fill=PCR))+
  geom_boxplot()+
  scale_fill_manual(values = c(cols$light_blue, cols$light_red))+
  stat_compare_means()+
  guides(fill="none")+
  labs(x="Covid", y="Error")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.12)))
ggsave(paste0(paths$out, "box_covid_pred_error.pdf"), height = 48, width = 63, units="mm", dpi = 300)
ggplot(meta_covid2, aes(x=PCR, y=Error_corrected, fill=PCR))+
  geom_boxplot()+
  scale_fill_manual(values = c(cols$light_blue, cols$light_red))+
  stat_compare_means()+
  guides(fill="none")+
  labs(x="Covid", y="Age-adjusted prediction error")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.12)))
ggsave(paste0(paths$out, "box_covid_pred_error_adj.pdf"), height = 48, width = 63, units="mm", dpi = 300)

sink(paste0(paths$out, "covid_tests.log"), type = "output")
t.test(Error~PCR, data=meta_covid2)
cat("=====================\n")
summary(lm(formula = Error~PCR, data=meta_covid2))
cat("=====================\n")
summary(lm(formula = Error~PCR+Age, data=meta_covid2))
closeAllConnections()

