library(matrixStats)
library(tidyverse)
library(data.table)
library(vegan)
library(patchwork)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggrepel)
library(gtools)
library(ggridges)

setwd("/scratch/fmorandi/ChromAcc-clock")

source("./scripts/utils.R")

##### LOAD DATA #####

paths = list()
paths$data = "./data/paper_data/"
paths$clocks_main = "./clocks/parallel/2023-02-03_14-38_tpm/"
# RNA vs ATAC
paths$clocks_atac_vs_rna = "./clocks/parallel/2023-02-04_15-10_atac_vs_rna/"
paths$clocks_rna_vs_atac = "./clocks/parallel/2023-02-04_15-10_rna_vs_atac/"
paths$clocks_multiomic = "./clocks/parallel/2023-02-27_12-09_multiomic/"
# Sex correction
paths$clocks_sex_crct = "./clocks/parallel/2023-02-20_17-32_crct_sex_vs_none/"
paths$clocks_sex_none = "./clocks/parallel/2023-02-20_17-31_crct_none_vs_sex/"
# Cell composition correction
paths$clocks_ccomp_only = "./clocks/parallel/2023-03-15_20-08_ccomp_only/"
paths$clocks_ccomp_and_atac = "./clocks/parallel/2023-03-15_20-10_crct_none_plus_ccomp/"
paths$clocks_ccomp_none = "./clocks/parallel/2023-03-15_20-08_crct_none_vs_ccomp/"
paths$clocks_ccomp_crct = "./clocks/parallel/2023-03-15_20-11_crct_ccomp_vs_none/"
paths$clocks_ccomp_age_crct = "./clocks/parallel/2023-03-16_09-23_crct_age/"
paths$clocks_ccomp_crct_no_age = "./clocks/parallel/2023-03-15_20-13_crct_ccomp_no_age/"
paths$clocks_ccomp_crct_in_ncv = "./clocks/parallel/2023-03-15_20-33_ccomp_in_ncv/"

enhancer_links = "peregrine"

# Load clean count data
atac = fread(paste0(paths$data, "data_atac_tpm.tsv"))
atac = data.frame(atac, row.names = 1)
rna = fread(paste0(paths$data, "data_rna_tmm.tsv"))
rna = data.frame(rna, row.names = 1)
atac_re = read.table(paste0(paths$data, "data_atac_re_cpm.tsv"))
atac_re_ocr = read.table(paste0(paths$data, "data_atac_re_ocr_perc.tsv"))

# Log transform
atac_log = log1p(atac)
rna_log = log1p(rna)

# Load metadata
meta = read.table(paste0(paths$data, "meta_final.tsv"))
rownames(meta) = meta$Subject
passingQC_atac = meta$Subject[which(meta$PassesQC_atac)]
passingQC_rna = meta$Subject[which(meta$PassesQC_rna)]

# Get labels
ages_atac = meta$Age
names(ages_atac) = meta$Subject
ages_rna = meta$Age
names(ages_rna) = meta$Subject

# Load qc table
qc_atac = read.table(paste0(paths$data, "qc_summary.tsv"), header=T)
qc_atac$SRR_atac = rownames(qc_atac)
qc_atac = merge(qc_atac, meta[,c("Subject", "SRR_atac", "Age")], by="SRR_atac")
rownames(qc_atac) = qc_atac$Subject

# Remove QC failed samples
ages_atac = ages_atac[passingQC_atac]
ages_rna = ages_rna[passingQC_rna]
atac = atac[passingQC_atac, ]
atac_log = atac_log[passingQC_atac, ]
rna = rna[passingQC_rna, ]
rna_log = rna_log[passingQC_rna, ]
atac_re = atac_re[passingQC_atac, ]
atac_re_ocr = atac_re_ocr[passingQC_atac, ]
qc_atac = qc_atac[passingQC_atac, ]

# Load peak_info and gene_info
peak_info = read.table(paste0(paths$data, "peak_info.tsv"), header=T)
colnames(peak_info) = c("desc", "chr", "start", "end", "strand", "len")
gene_info = read.table(paste0(paths$data, "gene_info.tsv"), header=T)

# Load anno
anno = read.table(paste0(paths$data, "peak_anno_", enhancer_links, ".txt"), header=T)

# Load ncv preditions
ncv_pred = list()
ncv_pred$main = read.table(paste0(paths$clocks_main, "preds.tsv"), row.names=1)
ncv_pred$atac_vs_rna = read.table(paste0(paths$clocks_atac_vs_rna, "preds.tsv"), row.names=1)
ncv_pred$rna_vs_atac = read.table(paste0(paths$clocks_rna_vs_atac, "preds.tsv"), row.names=1)
ncv_pred$multiomic = read.table(paste0(paths$clocks_multiomic, "preds.tsv"), row.names=1)
ncv_pred$sex_none = read.table(paste0(paths$clocks_sex_none, "preds.tsv"), row.names=1)
ncv_pred$sex_crct = read.table(paste0(paths$clocks_sex_crct, "preds.tsv"), row.names=1)
ncv_pred$ccomp_only = read.table(paste0(paths$clocks_ccomp_only, "preds.tsv"), row.names=1)
ncv_pred$ccomp_and_atac = read.table(paste0(paths$clocks_ccomp_and_atac, "preds.tsv"), row.names=1)
ncv_pred$ccomp_none = read.table(paste0(paths$clocks_ccomp_none, "preds.tsv"), row.names=1)
ncv_pred$ccomp_crct = read.table(paste0(paths$clocks_ccomp_crct, "preds.tsv"), row.names=1)
ncv_pred$ccomp_age_crct = read.table(paste0(paths$clocks_ccomp_age_crct, "preds.tsv"), row.names=1)
ncv_pred$ccomp_crct_no_age = read.table(paste0(paths$clocks_ccomp_crct_no_age, "preds.tsv"), row.names=1)
ncv_pred$ccomp_crct_in_ncv = read.table(paste0(paths$clocks_ccomp_crct_in_ncv, "preds.tsv"), row.names=1)
for (pred in names(ncv_pred)) {
  colnames(ncv_pred[[pred]]) = c("Age", "Preds", "Fold")
}

# Load ncv summaries
ncv_sum = list()
ncv_sum$main = read.table(paste0(paths$clocks_main, "summary.tsv"), row.names=1)
ncv_sum$atac_vs_rna = read.table(paste0(paths$clocks_atac_vs_rna, "summary.tsv"), row.names=1)
ncv_sum$rna_vs_atac = read.table(paste0(paths$clocks_rna_vs_atac, "summary.tsv"), row.names=1)
ncv_sum$multiomic = read.table(paste0(paths$clocks_multiomic, "summary.tsv"), row.names=1)
ncv_sum$sex_crct = read.table(paste0(paths$clocks_sex_crct, "summary.tsv"), row.names=1)
ncv_sum$sex_none = read.table(paste0(paths$clocks_sex_none, "summary.tsv"), row.names=1)
ncv_sum$ccomp_only = read.table(paste0(paths$clocks_ccomp_only, "summary.tsv"), row.names=1)
ncv_sum$ccomp_and_atac = read.table(paste0(paths$clocks_ccomp_and_atac, "summary.tsv"), row.names=1)
ncv_sum$ccomp_none = read.table(paste0(paths$clocks_ccomp_none, "summary.tsv"), row.names=1)
ncv_sum$ccomp_crct = read.table(paste0(paths$clocks_ccomp_crct, "summary.tsv"), row.names=1)
ncv_sum$ccomp_age_crct = read.table(paste0(paths$clocks_ccomp_age_crct, "summary.tsv"), row.names=1, fill = T)
ncv_sum$ccomp_crct_no_age = read.table(paste0(paths$clocks_ccomp_crct_no_age, "summary.tsv"), row.names=1)
ncv_sum$ccomp_crct_in_ncv = read.table(paste0(paths$clocks_ccomp_crct_in_ncv, "summary.tsv"), row.names=1)

# Load main clock ncv coefs
ncv_coef_main = read.table(paste0(paths$clocks_main, "/ncv_coefs.tsv"), row.names=1)
colnames(ncv_coef_main) = colnames(atac)

# Load ccomp + atac ncv coefs
ncv_coef_ccomp_and_atac = read.table(paste0(paths$clocks_ccomp_and_atac, "ncv_coefs.tsv"), row.names=1)
colnames(ncv_coef_ccomp_and_atac) = c(colnames(atac), colnames(meta)[16:23])

# Load main clock final coefs
final_coef = read.table(paste0(paths$clocks_main, "/final_coefs.tsv"), 
                        row.names=1, sep="\t", header = T)
final_clock_intercept = final_coef[nrow(final_coef), "coef"]
final_coef = final_coef[-c(nrow(final_coef)), ]

# Load multiomic clock final coefs
multiomic_coef = read.table(paste0(paths$clocks_multiomic, "/final_coefs.tsv"), row.names=1)
rownames(multiomic_coef) = c(colnames(atac), colnames(rna))
multiomic_coef$type = c(rep("ATAC", ncol(atac)), rep("RNA", ncol(rna)))
colnames(multiomic_coef) = c("coef", "type")

##### SETUP OUTPUT #####

dat = Sys.time()
dat = gsub(" ", "_", dat)
dat = gsub(":", "-", dat)

# Create output folders
paths$out = paste0(paths$data, "outputs_", dat, "/")
dir.create(paths$out)

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
cols$blue = "#0248A2"
cols$light_blue = "#95C1E3"
cols$dark_blue = "#062D67"
cols$red = "#8A0605"
cols$light_red = "#D25C5F"
cols$dark_red = "#550D0A"
cols$gray = "#DADADA"
cols$dark_gray = "#585858"
cols$yellow = "#EBB500" 
cols$light_yellow = "#FEE38A"
fold_cols = c("#B1DCFA","#0171BD","#062D67",
              "#631C17","#AA1C1A","#E15D66",
              "#EBB500","#FEE38A","#A3A3A3",
              "#777373","#000000")

# Others lying around
#0171BD #DADADA #AA1C19 (Blue, Gray, Red)
#AA1C19 = blue used in many points

# Redirect stdout to file
sink(paste0(paths$out, "main_script_res.log"), type = "output")
cat("Using", enhancer_links, "links\n")
cat("=====================\n")

##### SAMPLE COUNTS #####

# Number of ATAC samples before outlier removal
cat("ATAC:\tall samples:\ttotal:\t", sum(!is.na(meta$SRR_atac)), "\n")
cat("ATAC:\tall samples:\tmales:\t", sum(!is.na(meta$SRR_atac) & meta$Sex == "M"), "\n")
cat("ATAC:\tall samples:\tfemales:\t", sum(!is.na(meta$SRR_atac) & meta$Sex == "F"), "\n")

# Number of ATAC samples after outlier removal
cat("ATAC:\tfiltered samples:\ttotal:\t", sum(meta$PassesQC_atac), "\n")
cat("ATAC:\tfiltered samples:\tmales:\t", sum(meta$PassesQC_atac & meta$Sex == "M"), "\n")
cat("ATAC:\tfiltered samples:\tfemales:\t", sum(meta$PassesQC_atac & meta$Sex == "F"), "\n")

# Number of ATAC peaks
cat("ATAC:\tnumber of peaks:\t", ncol(atac), "\n")

# Number of RNA samples before outlier removal
cat("RNA:\tall samples:\ttotal:\t", sum(!is.na(meta$SRR_rna)), "\n")
cat("RNA:\tall samples:\tmales:\t", sum(!is.na(meta$SRR_rna) & meta$Sex == "M"), "\n")
cat("RNA:\tall samples:\tfemales:\t", sum(!is.na(meta$SRR_rna) & meta$Sex == "F"), "\n")

# Number of RNA samples after outlier removal
cat("RNA:\tfiltered samples:\ttotal:\t", sum(meta$PassesQC_rna), "\n")
cat("RNA:\tfiltered samples:\tmales:\t", sum(meta$PassesQC_rna & meta$Sex == "M"), "\n")
cat("RNA:\tfiltered samples:\tfemales:\t", sum(meta$PassesQC_rna & meta$Sex == "F"), "\n")

# Number of RNA genes
cat("RNA:\tnumber of genes:\t", ncol(rna), "\n")

# Number of samples passing outlier removal both in ATAC and RNA
cat("ATAC+RNA:\tfiltered samples:\ttotal:\t", sum(meta$PassesQC), "\n")
cat("ATAC+RNA:\tfiltered samples:\tmales:\t", sum(meta$PassesQC & meta$Sex == "M"), "\n")
cat("ATAC+RNA:\tfiltered samples:\tfemales:\t", sum(meta$PassesQC & meta$Sex == "F"), "\n")

cat("=====================\n")

##### ATAC AND RNA PCA #####

# ATAC
atac_pcs = prcomp(atac_log, center=TRUE, scale=FALSE)
expl_var = atac_pcs$sdev^2 / sum(atac_pcs$sdev^2)
atac_pcs = merge(atac_pcs$x, meta, by.x=0, by.y="Subject")
p1 = ggplot(atac_pcs, aes(x=PC1, y=PC2, color=Age)) +
  geom_point()+
  scale_color_gradient2(low=cols$blue, high=cols$red, midpoint=median(ages_atac))+
  labs(x = sprintf("PC1 [%.1f%%]", 100*expl_var[1]),
       y = sprintf("PC2 [%.1f%%]", 100*expl_var[2])) +
  ggtitle("ATAC-seq")
# Corr test and Adonis to test age effect
cor.test(atac_pcs$PC1, atac_pcs$Age)
cor.test(atac_pcs$PC2, atac_pcs$Age)
atac_dists = dist(atac_log)
adonis2(atac_dists ~ Age, data=meta[passingQC_atac, ]) # strata=meta[passingQC_atac, "Sex"]

# RNA
rna_pcs = prcomp(rna_log, center=TRUE, scale=FALSE)
expl_var = rna_pcs$sdev^2 / sum(rna_pcs$sdev^2)
rna_pcs = merge(rna_pcs$x, meta, by.x=0, by.y="Subject")
p2 = ggplot(rna_pcs, aes(x=PC1, y=PC2, color=Age)) +
  geom_point()+
  scale_color_gradient2(low=cols$blue, high=cols$red, midpoint=median(ages_atac))+
  labs(x = sprintf("PC1 [%.1f%%]", 100*expl_var[1]),
       y = sprintf("PC2 [%.1f%%]", 100*expl_var[2])) +
  ggtitle("RNA-seq", )
# Corr test and Adonis to test age effect
cor.test(rna_pcs$PC1, rna_pcs$Age)
cor.test(rna_pcs$PC2, rna_pcs$Age)
rna_dists = dist(rna_log)
adonis2(rna_dists ~ Age, data=meta[passingQC_rna, ])

(p1 | p2) + plot_layout(guides="collect")
ggsave(paste0(paths$out, "pca_atac_rna.pdf"), height = 72.4, width = w, units="mm", dpi = 300)

cat("=====================\n")

##### CELL COMPOSITION CORR PLOTS #####

ctypes = c("Monocytes", "Granulocytes", "Lymphocytes", "T_Cells", "CD4_T_Cells", "CD8_T_Cells", "B_Cells", "NK_Cells")
ps = list()
for (ctype in ctypes) {
  cell_lm_coef = coef(lm(as.formula(paste(ctype, "~ Age")), data = meta))
  p1 = ggplot(meta, aes_string(x="Age", y=ctype)) +
    geom_point(color=cols$azure)+
    geom_abline(intercept=cell_lm_coef[1], slope=cell_lm_coef[2],color="#AA1C19")+
    stat_cor(label.y = max(meta[ctype]*1.1, na.rm=T), cor.coef.name="r")+
    scale_y_continuous(expand = expansion(mult=c(0.05, 0.1)))
  # p1
  ps[[ctype]] = p1
}

wrap_plots(ps)
ggsave(paste0(paths$out, "corr_cells.pdf"), height = 122.8, width = w, units="mm", dpi = 300)

for (ctype in ctypes) {
  cell_cor_test = cor.test(meta$Age, meta[[ctype]], method="pearson")
  cat("CELL COMP CORR:\t", ctype, ":\tr =", round(cell_cor_test$estimate, digits = 4), 
      "\tp =", round(cell_cor_test$p.value, digits = 4), "\n")
}

cat("=====================\n")

##### GLOBAL DEREPRESSION  #####

lm_coef = coef(lm(frip ~ Age, data = qc_atac))
ggplot(qc_atac, aes(x=Age, y=frip))+
  geom_point(color="#0171BD")+
  geom_abline(intercept = lm_coef[1], slope=lm_coef[2], color="#AA1C19")+
  stat_cor(label.y = max(qc_atac$frip*1.1, na.rm=T), cor.coef.name="r")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.1)))+
  labs(y="Reads within OCRs (%)")
ggsave(paste0(paths$out, "corr_frip.pdf"), height = 45.3, width = 49.4, units="mm", dpi = 300)

atac_re$Repetitive = rowSums(atac_re[, 1:17]) / 1e4
atac_re$Age = ages_atac[rownames(atac_re)]
lm_coef = coef(lm(Repetitive ~ Age, data = atac_re))
ggplot(atac_re, aes(x=Age, y=Repetitive))+
  geom_point(color="#0171BD")+
  geom_abline(intercept = lm_coef[1], slope=lm_coef[2], color="#AA1C19")+
  stat_cor(label.y = max(atac_re$Repetitive*1.1, na.rm=T), cor.coef.name="r")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.1)))+
  labs(y="Reads over repetitive regions (%)")
ggsave(paste0(paths$out, "corr_repetitive.pdf"), height = 45.3, width = 49.4, units="mm", dpi = 300)

merge(meta[, c("Subject", "Age")], atac_re_ocr, by=0) %>%
  mutate(Other = 100-OCR-Repetitive) %>%
  mutate(Subject = fct_reorder(Subject, Age)) %>%
  dplyr::select(-Age, -Row.names) %>%
  pivot_longer(cols=-Subject) %>%
  mutate(name = factor(name, levels=c("OCR", "Other", "Repetitive"))) %>%
  ggplot(., aes(x=Subject, y=value, fill=name))+
    geom_bar(position="fill", stat="identity", width=1)+
    scale_fill_manual(values=c("#0171BD","#DADADA","#AA1C19"))+
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_x_discrete(expand = expansion(add = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    labs(x="Samples sorted by age", y="Reads (%)", fill="Genomic region")
ggsave(paste0(paths$out, "stack_ocr_repetitive.pdf"), height = 45.4, width = w, units="mm", dpi = 300)

# Test
frip_cor_test =  cor.test(qc_atac$Age, qc_atac$frip, method="pearson")
cat("FRIP CORR:\tr =", round(frip_cor_test$estimate, digits = 4), 
    "\tp =", round(frip_cor_test$p.value, digits = 4), "\n")
rep_cor_test =  cor.test(atac_re$Age, atac_re$Repetitive, method="pearson")
cat("REPETITIVE CORR:\tr =", round(rep_cor_test$estimate, digits = 4), 
    "\tp =", round(rep_cor_test$p.value, digits = 4), "\n")

atac_re %>%
  dplyr::select(-Age, -Repetitive) %>%
  cor_test_matrix(., ages_atac, type = "pearson") %>%
  mutate(p.adj = p.adjust(pvals, method="BH")) %>%
  mutate(Family = str_replace_all(rownames(.), "\\.", "?")) %>%
  ggplot(., aes(x=Family, y=1, fill=r, size=-log(p.adj)))+
    geom_point(pch=21)+
    geom_point(pch="*", aes(x=Family, alpha=p.adj<0.1), color="white")+
    scale_fill_gradient2(low=cols$blue, high=cols$red, limits = c(-0.4, 0.4))+
    scale_alpha_manual(values = c(0, 1)) +
    ylim(0, 2) +
    theme(axis.text.x=element_text(angle = 45, hjust=1),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          legend.position="top")+
    guides(alpha="none")+
    labs(x="Repetitive element family")
ggsave(paste0(paths$out, "corr_re_families.pdf"), height = 52.5, width = w, units="mm", dpi = 300)
    
cat("=====================\n")

##### ATAC SPEARMAN CORRELATION #####

sig_th_atac = 0.01
atac_age_cors = cor_test_matrix(atac_log, ages_atac, "spearman")

# FDR correction
atac_age_cors$padj = p.adjust(atac_age_cors$pvals, method = "fdr", n = nrow(atac_age_cors))
colnames(atac_age_cors) = c("corr","pval","padj")

# Save for integrative analysis and supplementary materials
peak_info = merge(peak_info, atac_age_cors, by=0, sort=F) %>%
  column_to_rownames("Row.names")
write.table(peak_info, paste0(paths$out, "peak_info_cors.tsv"), 
            sep = "\t", quote=F)

atac_age_cors_sig = atac_age_cors[atac_age_cors$padj < sig_th_atac,] 
cat("ATAC SPEAR:\tsignificant:\t", nrow(atac_age_cors_sig), "\n") # ?, 3037

atac_age_cors_pos = atac_age_cors_sig[atac_age_cors_sig$corr > 0,] 
cat("ATAC SPEAR:\tpositive:\t", nrow(atac_age_cors_pos), "\n") # 1113, 1124

atac_age_cors_neg = atac_age_cors_sig[atac_age_cors_sig$corr < 0,] 
cat("ATAC SPEAR:\tnegative:\t", nrow(atac_age_cors_neg), "\n") # 1891, 1913

atac_age_cors$sig = atac_age_cors$padj < sig_th_atac
atac_age_cors$sig_pos = atac_age_cors$sig & atac_age_cors$corr > 0
atac_age_cors$sig_neg = atac_age_cors$sig & atac_age_cors$corr < 0
atac_age_cors[, "sig_str"] = "Not sig."
atac_age_cors[atac_age_cors$sig_pos, "sig_str"] = "+"
atac_age_cors[atac_age_cors$sig_neg, "sig_str"] = "-"

# Funky bin choice so colored tails are straight
x1 = min(atac_age_cors$corr)
x2 = max(atac_age_cors[atac_age_cors$sig_neg, "corr"])
x3 = min(atac_age_cors[atac_age_cors$sig_pos, "corr"])
x4 = max(atac_age_cors$corr)
ncentral = ceiling(100 * (x3-x2) / (x4-x1))
binsize = (x3-x2) / ncentral
bins = seq(x1, x4, binsize)
closest_break = which.min(abs(bins - x2))
bins = bins - (bins - x2)[closest_break] * sign((bins - x2)[closest_break])
rm(x1, x2, x3, x4, ncentral, binsize, closest_break)

ggplot(atac_age_cors, aes(x=corr, fill=sig_str)) +
  geom_histogram(breaks=bins)+
  scale_fill_manual(values=c(cols$blue, cols$red, cols$gray))+
  guides(fill="none") +
  labs(x="r(ATAC-age)", y="Number of OCRs")
ggsave(paste0(paths$out, "hist_atac_spearman.pdf"), height = 51.8, width = 45.8, units="mm", dpi = 300)

cat("=====================\n")

##### SCATTER OF TOP CORRELATED OCRS #####

most_corr_pos = atac_age_cors %>%
  merge(., anno, by.x=0, by.y="peakID") %>%
  drop_na(gene_symbol) %>%
  dplyr::filter(gene_ensemblID %in% colnames(rna)) %>%
  top_n(5, corr) %>%
  arrange(-corr, desc(element_type)) %>%
  distinct(Row.names, .keep_all = T)

most_corr_neg = atac_age_cors %>%
  merge(., anno, by.x=0, by.y="peakID") %>%
  drop_na(gene_symbol) %>%
  dplyr::filter(gene_ensemblID %in% colnames(rna)) %>%
  top_n(5, -corr) %>%
  arrange(corr, desc(element_type)) %>%
  distinct(Row.names, .keep_all = T)

cat("Peaks with highest age correlations\n")
most_corr = rbind(most_corr_pos[1:2, ], 
                  most_corr_neg[1:2, ])
most_corr

ps = list()
for (i in 1:nrow(most_corr)) {
  titl = paste(most_corr$gene_symbol[i], most_corr$element_type[i])
  pk = most_corr$Row.names[i]
  tmp = data.frame(atac_log[pk], Age = ages_atac) # Rows match (I checked)
  # tmp$donor = rownames(tmp)
  lm_coef = coef(lm(as.formula(paste(pk, "~ Age")), data = tmp))
  ps[[i]] = ggplot(tmp, aes(x=Age, y=!!sym(pk)))+
    geom_point(color=cols$azure)+ #aes(color=donor %in% c("CR_036")) | color=cols$azure
    geom_abline(intercept=lm_coef[1], slope=lm_coef[2], color="#AA1C19")+
    stat_cor(label.y = max(tmp[pk]*1.1, na.rm=T), cor.coef.name="r")+
    scale_y_continuous(expand = expansion(mult=c(0.05, 0.1))) +
    labs(y="log(TPM)")+
    ggtitle(titl)+
    guides(color="none")
}
wrap_plots(ps, ncol=4)
ggsave(paste0(paths$out, "corr_top4.pdf"), height = 49, width = w, units="mm", dpi = 300)

cat("=====================\n")

##### PROMOTERS AND ENHANCERS AMONG CORRELATED OCRS #####

cat("element types: P&E, P_only, E_only, unannotated\n")

pos_peaks = element_type_lists(rownames(atac_age_cors_pos), anno)
pos_peaks_perc = 100 * unlist(pos_peaks[c("pro_and_enh_peaks_n", "pro_only_peaks_n", "enh_only_peaks_n", "unanno_peaks_n")]) / pos_peaks$all_peaks_n
cat("ATAC SPEAR:\tpositive:\telement_types:\t", pos_peaks_perc, "\n")

neg_peaks = element_type_lists(rownames(atac_age_cors_neg), anno)
neg_peaks_perc = 100 * unlist(neg_peaks[c("pro_and_enh_peaks_n", "pro_only_peaks_n", "enh_only_peaks_n", "unanno_peaks_n")]) / neg_peaks$all_peaks_n
cat("ATAC SPEAR:\tnegative:\telement_types:\t", neg_peaks_perc, "\n")

all_peaks = element_type_lists(rownames(atac_age_cors), anno)
all_peaks_perc = 100 * unlist(all_peaks[c("pro_and_enh_peaks_n", "pro_only_peaks_n", "enh_only_peaks_n", "unanno_peaks_n")]) / all_peaks$all_peaks_n
cat("ATAC SPEAR:\tall:\telement_types:\t", all_peaks_perc, "\n")

cat("=====================\n")

##### PROMOTERS AND ENHANCERS PIES AND STACKED BARS #####

labels = c("P&E","P","E","U")

data.frame(cbind(pos_peaks_perc, neg_peaks_perc, all_peaks_perc)) %>%
  mutate(type = labels) %>%
  pivot_longer(cols=-type)%>%
  mutate(name = factor(name, levels=c("neg_peaks_perc", "all_peaks_perc", "pos_peaks_perc"))) %>%
  ggplot(., aes(x=name, y=value,fill=name, alpha=type))+
    geom_bar(position="fill", stat="identity", color="black", linewidth = 0.2)+
    guides(fill="none") +
    scale_fill_manual(values = c(cols$dark_blue, "black", cols$dark_red)) +
    theme(legend.position="top", 
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    scale_x_discrete(labels=c("Closing", "All", "Opening"))+
    guides(fill="none", alpha="none")
ggsave(paste0(paths$out, "stack_element_type.pdf"), height = 26.5, width = 39.9, units="mm", dpi = 300)

##### PROMOTERS AND ENHANCERS ENRICHMENT TESTS AND PLOT #####

# Fisher's Exact and Chi Square Test

odds_ratios = c()
pvals = c()

# pos, promoters
cont_mat = make_contingency(pos_peaks$pro_peaks, pos_peaks$all_peaks, all_peaks$pro_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test: odds ratio = 0.472347, pval = 1.63e-29

# neg, promoters
cont_mat = make_contingency(neg_peaks$pro_peaks, neg_peaks$all_peaks, all_peaks$pro_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test: odds ratio= 0.6030838, pval= 4.55e-24

# pos, enhancers
cont_mat = make_contingency(pos_peaks$enh_peaks, pos_peaks$all_peaks, all_peaks$enh_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test: odds ratio = 2.62653, pval = 3.63e-48

# neg, enhancers
cont_mat = make_contingency(neg_peaks$enh_peaks, neg_peaks$all_peaks, all_peaks$enh_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test: odds ratio = 1.691909, pval = 4.75e-27

cat("enrichment: pos_pro, neg_pro, pos_enh, neg_enh\n")
cat("ATAC SPEAR:\tenrichment:\todds_ratios:\t", odds_ratios, "\n")
cat("ATAC SPEAR:\tenrichment:\tp_vals:\t", pvals, "\n")

data.frame(OR = odds_ratios, type = c("pro+", "pro-", "enh+", "enh-")) %>%
  mutate(type = factor(type, levels=c("pro-", "pro+", "enh-", "enh+"))) %>%
  ggplot(., aes(x=type, y=log(OR), fill=type))+
    geom_bar(stat="identity", color="black", linewidth = 0.2)+
    scale_fill_manual(values=c(cols$blue, cols$red, cols$light_blue, cols$light_red))+
    guides(fill="none")+
    labs(y="log(Odds ratio)")+
    theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
    scale_y_continuous(expand = expansion(mult=c(0.2, 0.2)))
ggsave(paste0(paths$out, "bar_element_type_enrich.pdf"), height = 19, width = 42.6, units="mm", dpi = 300)

cat("=====================\n")

##### RNA SPEARMAN CORRELATION #####

sig_th_rna = 0.01
rna_age_cors = cor_test_matrix(rna_log, ages_rna, "spearman")

# FDR correction
rna_age_cors$padj = p.adjust(rna_age_cors$pvals, method = "fdr", n = nrow(rna_age_cors))
colnames(rna_age_cors) = c("corr","pval","padj")

# Save for integrative analysis and supplementary materials
gene_info = merge(gene_info, rna_age_cors, by.x="Geneid", by.y=0, sort=F)
write.table(gene_info, paste0(paths$out, "gene_info_cors.tsv"), 
            sep = "\t", quote=F)

rna_age_cors_sig = rna_age_cors[rna_age_cors$padj < sig_th_rna,] 
cat("RNA SPEAR:\tsignificant:\t", nrow(rna_age_cors_sig), "\n") # 614

rna_age_cors_pos = rna_age_cors_sig[rna_age_cors_sig$corr > 0,] 
cat("RNA SPEAR:\tpositive:\t", nrow(rna_age_cors_pos), "\n") # 219

rna_age_cors_neg = rna_age_cors_sig[rna_age_cors_sig$corr < 0,] 
cat("RNA SPEAR:\tnegative:\t", nrow(rna_age_cors_neg), "\n") # 395

rna_age_cors$sig = rna_age_cors$padj < sig_th_rna
rna_age_cors$sig_pos = rna_age_cors$sig & rna_age_cors$corr > 0
rna_age_cors$sig_neg = rna_age_cors$sig & rna_age_cors$corr < 0
rna_age_cors[, "sig_str"] = "Not sig."
rna_age_cors[rna_age_cors$sig_pos, "sig_str"] = "+"
rna_age_cors[rna_age_cors$sig_neg, "sig_str"] = "-"

# Funky bin choice so colored tails are straight
x1 = min(rna_age_cors$corr)
x2 = max(rna_age_cors[rna_age_cors$sig_neg, "corr"])
x3 = min(rna_age_cors[rna_age_cors$sig_pos, "corr"])
x4 = max(rna_age_cors$corr)
ncentral = ceiling(100 * (x3-x2) / (x4-x1))
binsize = (x3-x2) / ncentral
bins = seq(x1, x4, binsize)
closest_break = which.min(abs(bins - x2))
bins = bins + (bins - x2)[closest_break] * sign((bins - x2)[closest_break])
rm(x1, x2, x3, x4, ncentral, binsize, closest_break)

ggplot(rna_age_cors, aes(x=corr, fill=sig_str)) +
  geom_histogram(breaks=bins)+
  scale_fill_manual(values=c(cols$blue, cols$red, cols$gray))+
  guides(fill="none") +
  labs(x="r(RNA-age)", y="Number of genes")
ggsave(paste0(paths$out, "hist_rna_spearman.pdf"), height = 0.3*h, width = 0.5*w, units="mm", dpi = 300)

cat("=====================/n")

##### GSEA: ATAC #####

gsea_links = anno[c("peakID", "gene_ensemblID")] %>%
  drop_na() %>%
  distinct() %>%
  dplyr::filter(gene_ensemblID %in% colnames(rna))
common_snames = intersect(rownames(atac_log), rownames(rna_log))
for (i in 1:nrow(gsea_links)) {
  pk = gsea_links[i, "peakID"]
  gn = gsea_links[i, "gene_ensemblID"]
  test = cor.test(rna_log[common_snames, gn], atac_log[common_snames, pk])
  gsea_links[i, "intercorr"] = test$estimate
  gsea_links[i, "pval"] = test$p.value
}
gsea_links = gsea_links %>% 
  group_by(gene_ensemblID) %>% 
  slice_max(order_by = intercorr)
# Convert ensembl to entrez, otherwise you sometimes get weirdness by gseGO
gsea_links$entrezID = mapIds(org.Hs.eg.db, keys = gsea_links$gene_ensemblID, 
                             keytype = "ENSEMBL", column = "ENTREZID") 
gsea_links = merge(gsea_links, atac_age_cors["corr"], by.x="peakID", by.y=0) %>%
  distinct(entrezID, .keep_all = T)

gene_list = gsea_links$corr
names(gene_list) = gsea_links$entrezID
gene_list = sort(gene_list, decreasing = TRUE)

set.seed(1337)
gse = gseGO(geneList=gene_list, 
            ont ="BP", 
            keyType = "ENTREZID",
            verbose = TRUE, 
            OrgDb = org.Hs.eg.db)
gsea_table = as.data.frame(gse)
# dotplot(gse, showCategory=6, split=".sign") +
#   facet_grid(.~.sign)+
#   theme(axis.text.y = element_text(size = 10))+
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
# gseaplot(gse, geneSetID="GO:0071216")

# Top 6 upregulated and top 6 downregulated by NES
gsea_table %>%
  arrange(NES) %>%
  slice_max(NES, n=6) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  mutate(log10_padj = -log10(p.adjust)) %>%
  ggplot(., aes(x=Description, y=NES, fill=log10_padj, size=1))+
    geom_point(pch=21)+
    scale_fill_gradient(high=cols$light_blue, low=cols$light_red)+
    coord_flip()+
    guides(size="none")
gsea_tails = rbind(slice_max(gsea_table, NES, n=6), slice_min(gsea_table, NES, n=6))
gsea_tails %>%
  mutate(Description = str_replace_all(Description, "tumor necrosis factor", "TNF")) %>%
  mutate(Description = str_replace_all(Description, "regulation", "reg.")) %>%
  mutate(Description = str_replace_all(Description, "heterochromatin", "heterochrom.")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(x=Description, y=NES, fill=NES>0))+
  geom_bar(stat="identity", linewidth=0.2, color="black", width = 0.8)+
  scale_fill_manual(values=c(cols$blue, cols$red))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 100))+
  guides(fill="none")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=6))+
  coord_flip()
# ggsave(paste0(paths$out, "bar_atac_gsea.pdf"), height = 51.6, width = 88.6, units="mm", dpi = 300)
ggsave(paste0(paths$out, "bar_atac_gsea_pere.pdf"), height = 51.6, width = w, units="mm", dpi = 300)

write.csv(gsea_table, paste0(paths$out, "gsea_atac.csv"))

##### GSEA: RNA #####

rna_age_cors$entrezID = mapIds(org.Hs.eg.db, keys = rownames(rna_age_cors), 
                               keytype = "ENSEMBL", column = "ENTREZID")
gene_list = rna_age_cors %>%
  drop_na() %>%
  distinct(entrezID, .keep_all = T) %>%
  arrange(-corr)%>%
  pull(corr, entrezID)

gc() # The second gseGO call sometimes causes some memory problems
set.seed(1337)
gse = gseGO(geneList=gene_list, 
            ont ="BP", 
            keyType = "ENTREZID",
            verbose = TRUE, 
            OrgDb = org.Hs.eg.db)
gsea_table = as.data.frame(gse)
# dotplot(gse, showCategory=20, split=".sign") + 
#   facet_grid(.~.sign)+ 
#   theme(axis.text.y = element_text(size = 10))+
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

gsea_tails = rbind(slice_max(gsea_table, NES, n=6), slice_min(gsea_table, NES, n=6))
gsea_tails %>%
  mutate(Description = str_replace_all(Description, "immunoglobulin", "Ig")) %>%
  mutate(Description = str_replace_all(Description, "mediated by", "by")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(x=Description, y=NES, fill=NES>0))+
  geom_bar(stat="identity", width=0.8, linewidth=0.2, color="black")+
  scale_fill_manual(values=c(cols$blue, cols$red))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 100))+
  guides(fill="none")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=6))+
  coord_flip()
ggsave(paste0(paths$out, "bar_rna_gsea.pdf"), height = 52.4, width = 88.6, units="mm", dpi = 300)

write.csv(gsea_table, paste0(paths$out, "gsea_rna.csv"))

##### ATAC-RNA INTEG: SIGNIFICANCE PLOT #####

# Annotation table without rows that have a missing gene link
anno_gene_nona = anno[!is.na(anno$gene_symbol),]
# Calculate -log(padj)
atac_age_cors$log.padj = -log(atac_age_cors$padj)
rna_age_cors$log.padj = -log(rna_age_cors$padj)
# Merge anno, atac_age_cors and rna_age_cors
atac_rna_pairs = merge(anno_gene_nona, atac_age_cors, by.x = "peakID", by.y = "row.names")
atac_rna_pairs = merge(atac_rna_pairs, rna_age_cors, by.x = "gene_ensemblID", by.y = "row.names", suffixes=c("_atac", "_rna"))

# Merged table only for atac-rna pairs where both are significant
atac_rna_pairs_sig = atac_rna_pairs %>%
  dplyr::filter(sig_atac & sig_rna) %>%
  mutate(color=interaction(element_type, ifelse(corr_atac > 0, "+", "-")))
ggplot(atac_rna_pairs_sig, aes(x=log.padj_atac, y=log.padj_rna, color=color))+
  geom_point()+
  geom_text_repel(aes(label=gene_symbol), max.overlaps=7, size=7*0.35)+
  scale_color_manual(values = c("#A9CEE9", "#0171BD","#DC7679", "#AA1C19"))+
  labs(x="-log(padj) of ATAC-Age correlation",
       y="-log(padj) of RNA-Age correlation") +
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.6)))+
  scale_x_continuous(expand = expansion(mult=c(0.05, 0.2)))+
  guides(color="none")
ggsave(paste0(paths$out, "scatter_atac_rna_sig.pdf"), height = 57.6, width = 86.9, units="mm", dpi = 300)

##### MAIN CLOCK NCV PERFORMANCE #####

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$main %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("ATAC clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("ATAC clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("ATAC clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("ATAC clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

# Annotation text for plot
anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$main %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_main.pdf"), height = 60.7, width = w/2, units="mm", dpi = 300)

cat("=====================\n")
  
##### MAIN CLOCK ELEMENT TYPES #####

clock_pos_peaks = rownames(final_coef)[final_coef$coef > 0]
clock_neg_peaks = rownames(final_coef)[final_coef$coef < 0]
clock_all_peaks = c(clock_pos_peaks, clock_neg_peaks)

cat("ATAC clock: final model selects", length(clock_all_peaks), "features,",
    length(clock_pos_peaks), "have a positive coef,",
    length(clock_neg_peaks), "have a negative coef\n")

cat("=====================\n")

clock_all_peaks = element_type_lists(clock_all_peaks, anno)
clock_all_peaks_perc = 100 * unlist(clock_all_peaks[c("pro_and_enh_peaks_n", "pro_only_peaks_n", "enh_only_peaks_n", "unanno_peaks_n")]) / clock_all_peaks$all_peaks_n

cat("element types: P&E, P_only, E_only, unannotated\n")
cat("ATAC clock:\telement_types:\t", clock_all_peaks_perc, "\n")

# Fisher's Exact and Chi Square Test

odds_ratios = c()
pvals = c()

# Promoters
cont_mat = make_contingency(clock_all_peaks$pro_peaks, clock_all_peaks$all_peaks, all_peaks$pro_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test Promoters: odds ratio = 0.7558744, pval = 0.02371

# Enhancers
cont_mat = make_contingency(clock_all_peaks$enh_peaks, clock_all_peaks$all_peaks, all_peaks$enh_peaks, all_peaks$all_peaks)
odds_ratios = append(odds_ratios, fisher.test(cont_mat)$estimate)
pvals = append(pvals, fisher.test(cont_mat)$p.value)
# Fisher's Exact Test: odds ratio = 1.114473, pval = 0.3788

cat("enrichment: promoters, enhancers (positive and negative coefs not separate)\n")
cat("ATAC clock:\tenrichment:\todds_ratios:\t", odds_ratios, "\n")
cat("ATAC clock:\tenrichment:\tp_vals:\t", pvals, "\n")

cat("=====================\n")

##### GSEA: CLOCK SITES #####

# No enrichment, skipping

##### TOP CLOCK SITES #####

# As opposed to atac_rna_pairs this keeps peaks that have no anno/link to no genes
pairs_all = merge(anno, atac_age_cors, by.x = "peakID", by.y = "row.names")
pairs_all = merge(pairs_all, rna_age_cors, by.x = "gene_ensemblID", by.y = "row.names", all.x = T,
                  suffixes = c("_atac", "_rna")) # <!>
pairs_all = merge(pairs_all, final_coef, by.x = "peakID", by.y = "row.names")

pairs_coef = pairs_all %>%
  dplyr::filter(!is.na(gene_symbol)) %>% # exclude peaks that have no linked gene (removes none if link to closest)
  # dplyr::filter(!is.na(corr_rna)) %>% # exclude peaks linked to non-expressed gene
  dplyr::filter(coef != 0) %>%
  arrange(desc(abs(coef))) %>%
  head(18) %>%
  dplyr::select(peakID, gene_symbol, element_type, corr_atac, padj_atac, corr_rna, padj_rna, coef) %>%
  unite("test", peakID:gene_symbol, remove = F) %>%
  mutate(dupe = duplicated(test))
pairs_coef[pairs_coef$dupe, "element_type"] = "Both"
pairs_coef = pairs_coef %>%
  arrange(peakID, gene_symbol, -dupe) %>%
  distinct_at(vars(-element_type, -dupe),.keep_all = T) %>%
  mutate(ord = rank(coef)) %>%
  arrange(ord) %>%
  mutate(ord=as.factor(ord))
pairs_coef[c("peakID", "gene_symbol", "element_type", "corr_atac", "padj_atac", "corr_rna", "padj_rna", "coef")]
p1 = pairs_coef %>%
  mutate(is_enhancer = element_type %in% c("Enhancer", "Both")) %>%
  mutate(is_promoter = element_type %in% c("Promoter", "Both")) %>%
  ggplot(.) +
    geom_point(aes(x=ord, y="Promoter", fill=is_promoter), pch=21, size=2)+
    geom_point(aes(x=ord, y="Enhancer", fill=is_enhancer), pch=21, size=2)+
    scale_fill_manual(values=c(cols$gray, "#000000"))+
    guides(fill="none", size="none")+
    coord_equal() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank())
p2 = pairs_coef %>%
  pivot_longer(cols=c("corr_atac", "corr_rna")) %>%
  ggplot(., aes(ord, y=name, fill=value))+
    geom_tile(aes(width=0.9, height=0.9))+
    scale_fill_gradient2(low=cols$blue, high=cols$red) +
    coord_equal() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank())
p3 = ggplot(pairs_coef, aes(x=ord, y="Coef", fill=coef))+
  geom_tile(aes(width=0.9, height=0.9))+
  scale_fill_gradient2(low=cols$blue, high=cols$red) +
  coord_equal() +
  theme(axis.text.x=element_text(angle = 45, hjust=1),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  scale_x_discrete(labels=pairs_coef$gene_symbol)
p1/p2/p3
ggsave(paste0(paths$out, "top_clock_sites.pdf"), height = 57.9, width = 93.9, units="mm", dpi = 300)

##### CLOCK COMPAPRISON: ATAC VS RNA #####

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$atac_vs_rna %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("ATAC clock vs RNA: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("ATAC clock vs RNA: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("ATAC clock vs RNA: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("ATAC clock vs RNA: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

# Annotation text for plot
anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$atac_vs_rna %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_atac_vs_rna.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$rna_vs_atac %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("RNA clock vs ATAC: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("RNA clock vs ATAC: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("RNA clock vs ATAC: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("RNA clock vs ATAC: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

# Annotation text for plot
anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$rna_vs_atac %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_rna_vs_atac.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$multiomic %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Multiomic clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Multiomic clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Multiomic clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Multiomic clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

# Annotation text for plot
anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$multiomic %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_multiomic.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Correlation of clock predictions
aVr_agreement = cor.test(ncv_pred$atac_vs_rna$Preds, ncv_pred$rna_vs_atac$Preds)
cat("Corr of ATAC and RNA clock pred: r =", aVr_agreement$estimate, ", p =", aVr_agreement$p.value, "\n")

# Performance comparison boxplot
ncv_sum$atac_vs_rna$type = "ATAC"
ncv_sum$rna_vs_atac$type = "RNA"
ncv_sum$multiomic$type = "Multiomic"
atac_rna_perf = rbind(ncv_sum$atac_vs_rna, ncv_sum$rna_vs_atac, ncv_sum$multiomic) %>%
  dplyr::select(-log10.alpha., -l1_ratio, -non_zero) %>%
  pivot_longer(cols=-type) %>%
  mutate(name = factor(name, levels=c("RMSE", "MAE", "r")))
ggplot(atac_rna_perf, aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c(cols$light_yellow, cols$yellow, cols$gray))+
  facet_wrap(~name, scales="free")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(fill = "")
ggsave(paste0(paths$out, "box_atac_vs_rna_vs_multiomic.pdf"), height = 58.1, width = 96.6, units="mm", dpi = 300)

compare_means(value~type, data=atac_rna_perf, group.by="name", method="t.test")

multiomic_coef %>%
  dplyr::filter(coef != 0) %>%
  ggplot(., aes(x=type, y = abs(coef), fill=type))+
    geom_boxplot()+
    scale_fill_manual(values=c(cols$red, cols$blue))+
    stat_compare_means(method="t.test")
ggsave(paste0(paths$out, "box_multiomic_coefs.pdf"), height = 61.2, width = 57.8, units="mm", dpi = 300)

cat("=====================\n")

fg = multiomic_coef %>%
  dplyr::filter(coef != 0) %>%
  group_by(type) %>%
  summarize(in_clock = n()) %>%
  column_to_rownames(var="type")
bg = multiomic_coef %>%
  dplyr::filter(coef == 0) %>%
  group_by(type) %>%
  summarize(in_bg = n()) %>%
  column_to_rownames(var="type")
cont_table = cbind(fg, bg)
cont_table
fisher.test(cont_table)

##### CLOCK COMPAPRISON: SEX CORRECTION #####

ncv_pred$main %>% 
  merge(., meta, by=0) %>%
  mutate(error = Preds-Age.x) %>%
  ggplot(., aes(x=Sex, y=error, fill=Sex))+
  geom_boxplot()+
  scale_fill_manual(values=c(cols$light_red, cols$light_blue))+
  stat_compare_means(method="t.test")
ggsave(paste0(paths$out, "box_sex_main_clock_error.pdf"), height = 50.5, width = 47.8, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$sex_crct %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Sex corrected clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Sex corrected clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Sex corrected clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Sex corrected clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$sex_crct %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_sex_crct.pdf"), height = 50.4, width = 58.4, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$sex_none %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Sex uncorrected clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Sex uncorrected clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Sex uncorrected clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Sex uncorrected clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$sex_none %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_sex_none.pdf"), height = 50.4, width = 58.4, units="mm", dpi = 300)

cat("=====================\n")

# Performance comparison boxplot
ncv_sum$sex_crct$type = "sex_crct"
ncv_sum$sex_none$type = "sex_none"
sex_perf = rbind(ncv_sum$sex_crct, ncv_sum$sex_none) %>%
  dplyr::select(-log10.alpha., -l1_ratio, -non_zero) %>%
  pivot_longer(cols=-type) %>%
  mutate(name = factor(name, levels=c("RMSE", "MAE", "r")))
ggplot(sex_perf, aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c(cols$light_yellow, cols$yellow))+
  facet_wrap(~name, scales="free")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(fill = "")
ggsave(paste0(paths$out, "box_sex_crct_vs_none.pdf"), height = 58.1, width = 96.6, units="mm", dpi = 300)

compare_means(value~type, data=sex_perf, group.by="name", method="t.test")

cat("=====================\n")

##### CLOCK COMPAPRISON: CELL COMPOSITION #####

###### Are clock sites more correlated with age or cell composition? ######
atac_log_clock_sites = atac_log[, rownames(final_coef)[final_coef$coef != 0]]
ccomp_and_age = meta %>%
  dplyr::select(Age, Monocytes:CD8_T_Cells) %>%
  drop_na()

clock_site_cors = data.frame(
  cor(atac_log_clock_sites, ccomp_and_age[rownames(atac_log_clock_sites),],
      use="pairwise.complete.obs")) %>%
  rownames_to_column(var="peakID") %>%
  pivot_longer(cols = -peakID, names_to = "Variable", values_to = "r") %>%
  merge(., final_coef["coef"], by.x="peakID", by.y=0) %>%
  mutate(peakID = fct_reorder(peakID, coef)) %>%
  arrange(peakID)
max_cors = clock_site_cors %>%
  group_by(peakID) %>%
  slice_max(abs(r))
cat("Variable best correlating for each clock site:/n")
table(max_cors$Variable)
  
# Cors of all clock sites with age and ccomp
ggplot(clock_site_cors, aes(x=peakID, y=Variable, fill=r))+
  geom_tile()+
  scale_fill_gradient2(low=cols$blue, mid="white", high=cols$red)+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank())+
  scale_x_discrete(expand = expansion(add = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(0, 0)))
ggsave(paste0(paths$out, "heatmap_clock_site_cors.pdf"), height = 50.4, width = w, units="mm", dpi = 300)

# Density of cors of all clock sites with age and ccomp
ggplot(clock_site_cors, aes(x=r, y=Variable))+
  geom_density_ridges()
ggsave(paste0(paths$out, "ridges_clock_site_cors.pdf"), height = 50.4, width = 0.5*w, units="mm", dpi = 300)

cat("=====================\n")

######  Does cell composition alone allow for age prediction? ###### 
# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$ccomp_only %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Cell composition clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Cell composition clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Cell composition clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Cell composition clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$ccomp_only %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_ccomp_only.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

######  Does the clock ever pick cell composition features if provided besides atac? ###### 

cat("Number of ccomp features picked at least once:\n")
data.frame(t(ncv_coef_ccomp_and_atac)) %>%
  rownames_to_column(var="var") %>% 
  pivot_longer(cols = -var) %>%
  dplyr::filter(value!=0) %>%
  group_by(var) %>%
  summarize(n()) %>%
  dplyr::filter(!grepl("peak", var)) %>%
  nrow()

cat("=====================\n")

######  Compare uncorrected, corrected for ccomp, corrected for age ######
# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$ccomp_none %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Cell composition uncorrected clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Cell composition uncorrected clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Cell composition uncorrected clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Cell composition uncorrected clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$ccomp_none %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_ccomp_none.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$ccomp_crct %>%
  dplyr::select(non_zero:r) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Cell composition corrected clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Cell composition corrected clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Cell composition corrected clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Cell composition corrected clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$ccomp_crct %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_ccomp_crct.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Mean and sd of performance metrics
perf_mean_sd = ncv_sum$ccomp_age_crct %>%
  dplyr::select(non_zero:r) %>%
  mutate(r = replace_na(r, 0)) %>%
  rowid_to_column("Fold") %>%
  pivot_longer(cols=-Fold) %>%
  group_by(name) %>%
  summarize(mean = mean(value), sd = sd(value)) %>%
  column_to_rownames("name")

# Print to log
cat("Age corrected clock: n features =", round(perf_mean_sd["non_zero", "mean"], 0), "±", 
    round(perf_mean_sd["non_zero", "sd"], 0), "\n")
cat("Age corrected clock: RMSE =", round(perf_mean_sd["RMSE", "mean"], 2), "±", 
    round(perf_mean_sd["RMSE", "sd"], 2), "\n")
cat("Age corrected clock: MAE =", round(perf_mean_sd["MAE", "mean"], 2), "±", 
    round(perf_mean_sd["MAE", "sd"], 2), "\n")
cat("Age corrected clock: r =", round(perf_mean_sd["r", "mean"], 2), "±", 
    round(perf_mean_sd["r", "sd"], 2), "\n")

anno_text = sprintf("RMSE = %.2f ± %.2f\nMAE = %.2f ± %.2f\nr = %.2f ± %.2f", 
                    perf_mean_sd["RMSE", "mean"], perf_mean_sd["RMSE", "sd"],
                    perf_mean_sd["MAE", "mean"], perf_mean_sd["MAE", "sd"],
                    perf_mean_sd["r", "mean"], perf_mean_sd["r", "sd"])

# Performance scatter plot
ncv_pred$ccomp_age_crct %>%
  mutate(Fold = as.factor(Fold)) %>%
  ggplot(., aes(x=Age, y=Preds, color=Fold))+
  geom_point()+
  geom_abline(slope=1, color=cols$dark_red)+
  scale_color_manual(values=fold_cols)+
  labs(y="Predicted Age")+
  scale_x_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  scale_y_continuous(breaks = seq(20, 73, by = 10), limits = c(19, 74))+
  annotate("text", x = 19, y=74, label=anno_text, size=7*0.35, vjust=1, hjust=0)+
  guides(color="none")
ggsave(paste0(paths$out, "corr_clock_age_crct.pdf"), height = 50.4, width = 52, units="mm", dpi = 300)

cat("=====================\n")

# Performance comparison boxplot
ncv_sum$ccomp_only$type = "ccomp_only"
ncv_sum$ccomp_none$type = "ccomp_none"
ncv_sum$ccomp_crct$type = "ccomp_crct"
ncv_sum$ccomp_age_crct$type = "age_crct"

ccomp_perf = rbind(ncv_sum$ccomp_only, ncv_sum$ccomp_none, 
                   ncv_sum$ccomp_crct, ncv_sum$ccomp_age_crct) %>%
  dplyr::select(-log10.alpha., -l1_ratio, -non_zero) %>%
  pivot_longer(cols=-type) %>%
  mutate(name = factor(name, levels=c("RMSE", "MAE", "r"))) %>%
  mutate(type = factor(type, levels=c("ccomp_only", "ccomp_none", "ccomp_crct", "age_crct")))
ggplot(ccomp_perf, aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c(cols$light_yellow, cols$gray, cols$yellow, cols$dark_gray))+
  facet_wrap(~name, scales="free")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(fill = "")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.2)))
ggsave(paste0(paths$out, "box_ccomp_crct1.pdf"), height = 58.1, width = w, units="mm", dpi = 300)

compare_means(value~type, data=ccomp_perf, group.by="name", method="t.test")

cat("=====================\n")

###### Can we correct for cell composition without age info? ######

# Performance comparison boxplot
ncv_sum$ccomp_none$type = "ccomp_none"
ncv_sum$ccomp_crct$type = "ccomp_crct"
ncv_sum$ccomp_crct_no_age$type = "ccom_crct_no_age"
ncv_sum$ccomp_crct_in_ncv$type = "ccomp_crct_in_ncv"

ccomp_perf = rbind(ncv_sum$ccomp_none, ncv_sum$ccomp_crct, 
                   ncv_sum$ccomp_crct_no_age, ncv_sum$ccomp_crct_in_ncv) %>%
  dplyr::select(-log10.alpha., -l1_ratio, -non_zero) %>%
  pivot_longer(cols=-type) %>%
  mutate(name = factor(name, levels=c("RMSE", "MAE", "r"))) %>%
  mutate(type = factor(type, levels=c("ccomp_none", "ccom_crct_no_age", "ccomp_crct_in_ncv", "ccomp_crct")))
ggplot(ccomp_perf, aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c(cols$light_yellow, cols$gray, cols$yellow, cols$dark_gray))+
  facet_wrap(~name, scales="free")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(fill = "")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.2)))
ggsave(paste0(paths$out, "box_ccomp_crct2.pdf"), height = 58.1, width = w, units="mm", dpi = 300)

compare_means(value~type, data=ccomp_perf, group.by="name", method="t.test")

closeAllConnections()
