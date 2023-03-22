library(tidyverse)
library(limma)
library(edgeR)
library(sp)
library(patchwork)
library(scales)

setwd("/scratch/fmorandi/ChromAcc-clock")

paths = list()
paths$pipeline ="./data/pipeline_outputs/atac/"
paths$meta="./data/meta.tsv"
paths$out="./data/paper_data/"

# Options
th_counts = 1.1e7
th_frip = 0.18
w = 174 # mm
h = 230

dir.create(paste0(paths$out, "figures"), showWarnings = F)

##### LOAD DATA #####

counts = read.table(paste0(paths$pipeline, "07_peaksets_and_tables/counts.tsv"), header=T)
peak_info = counts[, 1:6]
counts = counts[, -c(1:6)]
colnames(counts) = str_extract(colnames(counts), "(SRR[0-9]+)", group=1)
rownames(peak_info) = paste0("peak", 1:nrow(peak_info))
rownames(counts) = rownames(peak_info)
counts = data.frame(t(counts))

fc_report = read.table(paste0(paths$pipeline, "07_peaksets_and_tables/counts.tsv.summary"),
                       header=T, row.names=1)
colnames(fc_report) = str_extract_all(colnames(fc_report), "(SRR[0-9]+)")

meta = read.table(paste0(paths$meta), header=T)

##### EXTRACT QC INFO #####

logs = dir(paste0(paths$pipeline, "A_mapping_logs_old"), full.names=T)

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
all(colnames(fc_report) == rownames(qc))
all(colSums(fc_report)/2 == qc$clean)
# Get frip
fc_report = data.frame(t(fc_report))
qc$cut_sites = fc_report$Assigned + fc_report$Unassigned_NoFeatures
qc$cut_sites_in_peak = fc_report$Assigned
qc$frip = qc$cut_sites_in_peak / qc$cut_sites

##### CONVERT NAMES #####

dict = data.frame(Subject = meta$Subject, SRR = meta$SRR_atac)
dict = drop_na(dict)
rownames(dict) = dict$SRR

rownames(counts) = dict[rownames(counts), "Subject"]

##### MARK LOW QUALITY SAMPLES #####

ggplot(qc, aes(x=clean, y=frip))+
  geom_point()+
  geom_vline(xintercept=th_counts, color="red")+
  geom_hline(yintercept=th_frip, color="red")+
  labs(x="N of filtered pairs", y="FRIP") +
  theme_light()
ggsave(paste0(paths$out, "figures/qc_atac_thresholds.pdf"), height = 0.3*h, width = 0.5*w, units="mm", dpi = 300)

bad_samples = rownames(qc)[qc$clean < th_counts | qc$frip < th_frip]

meta$LowQ_atac = meta$SRR_atac %in% bad_samples

##### TPM NORMALIZATION #####

tpm = 1000 * sweep(counts, 2, peak_info$Length, "/")
tpm = 1e6 * sweep(tpm, 1, rowSums(tpm), "/")

##### REMOVE OUTLIERS #####

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
ggsave(paste0(paths$out, "figures/qc_atac_ellipse.pdf"), height = 0.4*h, width = 1*w, units="mm", dpi = 300)

meta$Outlier_atac = meta$Subject %in% pca$Row.names[pca$outlier]
meta$PassesQC_atac = !is.na(meta$SRR_atac) & !meta$LowQ_atac & !meta$Outlier_atac

##### SEX AND CELL COMPOSITION CORRECTION #####

meta2 = meta
meta2 = merge(meta2, qc, by.x="SRR_atac", by.y=0)
meta2 = meta2[meta2$PassesQC_atac, ]
rownames(meta2) = meta2$Subject

dge = DGEList(counts=t(counts[rownames(meta2), ]), samples = meta2)
dge = calcNormFactors(dge, method="TMM")

# No correction
design = model.matrix(~Age, data=dge$samples)
v = voom(dge, design, plot=T)
voom_crct_none = data.frame(t(as.matrix(v)))

pca = prcomp(voom_crct_none)
pca = pca$x[, 1:10]
all(rownames(pca) == rownames(meta2))
vars = c("Age", "Sex", "Monocytes", "Granulocytes", "Lymphocytes",
        "B_Cells", "NK_Cells", "T_Cells", "CD4_T_Cells", "CD8_T_Cells",
        "mitochondrial", "duplicates", "clean", "frip")
meta2 = meta2[vars]
meta2$Sex = as.integer(as.factor(meta2$Sex))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m1 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("No correction")
m1

# Sex correction
design = model.matrix(~Age, data=dge$samples)
voom_crct_sex = removeBatchEffect(v, batch=dge$samples$Sex, design=design)
voom_crct_sex = data.frame(t(voom_crct_sex))

pca = prcomp(voom_crct_sex)
pca = pca$x[, 1:10]
all(rownames(pca) == rownames(meta2))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m2 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Correction for sex")
m2

# Cell composition correction knowing age
dge = dge[, which(!is.na(dge$samples$Monocytes))]
v = v[, colnames(dge)]
design = model.matrix(~Age, data=dge$samples)
voom_crct_ccomp = removeBatchEffect(v, covariates=dge$samples[19:26], design=design)
voom_crct_ccomp = data.frame(t(voom_crct_ccomp))

pca = prcomp(voom_crct_ccomp)
pca = pca$x[, 1:10]
meta2 = dge$samples[vars]
meta2$Sex = as.integer(as.factor(meta2$Sex))
all(rownames(pca) == rownames(meta2))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m3 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Correction for ccomp, knowing age")
m3

# Sex and cell composition correction knowing age
voom_crct_both = removeBatchEffect(v, batch = dge$samples$Sex, covariates=dge$samples[19:26], design=design)
voom_crct_both = data.frame(t(voom_crct_both))

pca = prcomp(voom_crct_both)
pca = pca$x[, 1:10]
meta2 = dge$samples[vars]
meta2$Sex = as.integer(as.factor(meta2$Sex))
all(rownames(pca) == rownames(meta2))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m4 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Correction for ccomp and sex, knowing age")
m4

# Cell composition correction with no age info
voom_crct_ccomp_no_age = removeBatchEffect(v, covariates=dge$samples[19:26])
voom_crct_ccomp_no_age = data.frame(t(voom_crct_ccomp_no_age))

pca = prcomp(voom_crct_ccomp_no_age)
pca = pca$x[, 1:10]
meta2 = dge$samples[vars]
meta2$Sex = as.integer(as.factor(meta2$Sex))
all(rownames(pca) == rownames(meta2))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m5 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Correction for ccomp, not knowing age")
m5

# Cell composition correction with no age info
design = model.matrix(~Monocytes+Granulocytes+Lymphocytes+B_Cells+
                        NK_Cells+T_Cells+CD4_T_Cells+CD8_T_Cells, data=dge$samples)
voom_crct_age = removeBatchEffect(v, covariates=dge$samples$Age, design=design)
voom_crct_age = data.frame(t(voom_crct_age))

pca = prcomp(voom_crct_age)
pca = pca$x[, 1:10]
meta2 = dge$samples[vars]
meta2$Sex = as.integer(as.factor(meta2$Sex))
all(rownames(pca) == rownames(meta2))
pca_cors = cor(pca, meta2, use="pairwise.complete.obs")
m6 = data.frame(pca_cors) %>%
  mutate(PC = rownames(.)) %>%
  pivot_longer(cols=-PC, names_to = "Variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels=paste0("PC", 1:10))) %>%
  mutate(Variable = factor(Variable, levels=vars)) %>%
  ggplot(., aes(x=Variable, y=PC, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Correction for age, leaving ccomp")
m6

ps = align_patches(m1, m2, m3, m4, m5, m6)
pdf(paste0(paths$out, "figures/qc_atac_correction.pdf"))
for (p in ps) {
  plot(p)
}
dev.off()

##### SUMMARIZE CLEANUP #####

p1 = ggplot(meta, aes(x=Age, fill=!PassesQC_atac))+
  geom_histogram(breaks=0.5+seq(20, 74))+
  labs(fill="Failing QC") +
  theme_light()
p2 = ggplot(meta, aes(x=Sex, fill=!PassesQC_atac))+
  geom_bar(position="fill")+
  labs(fill="Failing QC") +
  theme_light()
(p1 | p2) + plot_layout(widths = c(4, 1), guides = 'collect')
ggsave(paste0(paths$out, "figures/qc_atac_distrs.pdf"), height = 0.3*h, width = w, units="mm", dpi = 300)

##### SAVE OUTPUTS #####

writeLines(meta[!meta$PassesQC_atac, "Subject"], paste0(paths$out, "qc_failed_atac.txt"))

peak_bed = peak_info[, c("Chr", "Start", "End")]
peak_bed$Name = rownames(peak_info)

write.table(tpm, file=paste0(paths$out, "data_atac_tpm.tsv"), quote=F, sep="\t")
write.table(voom_crct_none, file=paste0(paths$out, "data_atac_crct_none.tsv"), quote=F, sep="\t")
write.table(voom_crct_sex, file=paste0(paths$out, "data_atac_crct_sex.tsv"), quote=F, sep="\t")
write.table(voom_crct_ccomp, file=paste0(paths$out, "data_atac_crct_ccomp.tsv"), quote=F, sep="\t")
write.table(voom_crct_both, file=paste0(paths$out, "data_atac_crct_both.tsv"), quote=F, sep="\t")
write.table(voom_crct_ccomp_no_age, file=paste0(paths$out, "data_atac_crct_ccomp_no_age.tsv"), quote=F, sep="\t")
write.table(voom_crct_age, file=paste0(paths$out, "data_atac_crct_age.tsv"), quote=F, sep="\t")

write.table(peak_info, file=paste0(paths$out, "peak_info.tsv"), quote=F, sep="\t")
write.table(peak_bed, file=paste0(paths$out, "peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(meta, file=paste0(paths$out, "meta_atac.tsv"), quote=F, sep="\t")
write.table(qc, file=paste0(paths$out, "qc_summary.tsv"), quote=F, sep="\t")

