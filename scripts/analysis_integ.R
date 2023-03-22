library(matrixStats)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("/scratch/fmorandi/ChromAcc-clock")

source("./scripts/utils.R")

##### LOAD DATA #####

paths = list()
paths$data = "./data/paper_data/"
paths$clocks_main = "./clocks/parallel/2023-02-03_14-38_tpm/"
paths$met = "./data/external_resources/methylation/"
paths$out = "./data/paper_data/outputs_2023-03-20_15-27-28_blood/"

enhancer_links = "blood"

# Load metadata
meta = read.table(paste0(paths$data, "meta_final.tsv"))
rownames(meta) = meta$Subject

# Load peak_info and gene_info
peak_info = read.table(paste0(paths$out, "peak_info_cors.tsv"), header=T)
gene_info = read.table(paste0(paths$out, "gene_info_cors.tsv"), header=T)

# Load anno
anno = read.table(paste0(paths$data, "peak_anno_", enhancer_links, ".txt"), header=T)

# Load main clock final coefs
final_coef = fread(paste0(paths$clocks_main, "final_coefs.tsv"), header = T)
final_coef = data.frame(final_coef, row.names = 1)
final_coef = final_coef[final_coef$coef != 0, ]
final_coef = final_coef[-which(rownames(final_coef) == "intercept"), ]

clock_pos_peaks = rownames(final_coef)[final_coef$coef > 0]
clock_neg_peaks = rownames(final_coef)[final_coef$coef < 0]
clock_all_peaks = c(clock_pos_peaks, clock_neg_peaks)

# Load Hannum methylation-age correlation table (pre-computed as it takes long)
hannum_age_cors = read.table(paste0(paths$data, "hannum_age_cors.tsv"), header = 1)

# Load cpg-peak mapping
# Note: there are some probe ids that are in cpg_peak_pairs but not in met_age_cors
# cpg_peak_pairs comes from the illumina list of probes, intersected to our peaks
# probably some were filtered by hannum et al. before uploading the table
cpg_peak_pairs = read.table(paste0(paths$data, "peak_cpg_map.tsv"), header = 1)

# Load Hannum and Horvath coefs
hannum_coef = read.csv(paste0(paths$met, "hannum_coef.csv"), header=TRUE, row.names = 1)
horvath_coef = read.csv(paste0(paths$met, "horvath_coef.csv"), header=TRUE, row.names = 1)

##### SETUP OUTPUT #####

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

# Redirect stdout to file
sink(paste0(paths$out, "integ_script_res.log"), type = "output")
cat("Using", enhancer_links, "links\n")
cat("=====================\n")

##### ATAC-MET INTEG: NUMBERS #####

el_type = anno %>%
  dplyr::group_by(peakID) %>%
  summarize(element_type = paste(unique(element_type), collapse = "-")) %>%
  mutate(element_type = fct_recode(element_type, "Both" = "Enhancer-Promoter", 
                                   "Both" = "Promoter-Enhancer", "Unanno" = "NA")) %>%
  merge(cpg_peak_pairs, by.x="peakID", by.y = "peak_id") %>%
  filter(probe_id %in% rownames(hannum_age_cors))

ncpg_pro = sum(el_type$element_type == "Promoter")
ncpg_enh = sum(el_type$element_type == "Enhancer")
ncpg_both = sum(el_type$element_type == "Both")
ncpg_unanno = sum(el_type$element_type == "Unanno")

cat("MET numbers:\t", nrow(cpg_peak_pairs), "cpgs out of", nrow(hannum_age_cors), "are inside our OCRs\n")
cat("MET numbers:\t", ncpg_pro, "lie in promoters\n")
cat("MET numbers:\t", ncpg_enh, "lie in enhancers\n")
cat("MET numbers:\t", ncpg_both, "lie in promoter-enhancers\n")
cat("MET numbers:\t", ncpg_unanno, "lie in unannotated OCRs\n")

cat("=====================\n")

##### HEATMAPS: RNA VS ATAC VS MET #####

triads = merge(anno, peak_info[,c("corr", "pval", "padj")], by.x = "peakID", by.y = "row.names")
# I want to keep all peaks even if they dont have gene links because they might contain CpGs
triads = merge(triads, gene_info[,c("Geneid", "corr", "pval", "padj")], by.x = "gene_ensemblID", by.y = "Geneid", 
               suffixes = c("_atac", "_rna"), all.x=T)
# I want to keep all peaks even if they dont contain probes because they might link to genes
triads = merge(triads, cpg_peak_pairs, by.x = "peakID", by.y = "peak_id", all.x=T)
triads = hannum_age_cors %>%
  dplyr::rename("corr_met"="corr", "pval_met" = "pval", "padj_met" = "padj") %>%
  merge(triads, ., by.x = "probe_id", by.y=0, all.x = T)
triads = triads %>%
  mutate(sig_rna = interaction(padj_rna < 0.01, sign(corr_rna)), .after = padj_rna) %>%
  mutate(sig_rna = fct_recode(sig_rna, Increasing = "TRUE.1",  Decreasing = "TRUE.-1", 
                              Constant = "FALSE.1", Constant = "FALSE.-1")) %>%
  mutate(sig_atac = interaction(padj_atac < 0.01, sign(corr_atac)), .after = padj_atac) %>%
  mutate(sig_atac = fct_recode(sig_atac, Increasing = "TRUE.1",  Decreasing = "TRUE.-1", 
                               Constant = "FALSE.1", Constant = "FALSE.-1", Constant = "FALSE.0")) %>%
  mutate(sig_met = interaction(padj_met < 0.01, sign(corr_met)), .after = padj_met) %>%
  mutate(sig_met = fct_recode(sig_met, Increasing = "TRUE.1",  Decreasing = "TRUE.-1", 
                              Constant = "FALSE.1", Constant = "FALSE.-1"))

# Promoters
triad_cors = triads %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("Promoters")
ggsave(paste0(paths$out, "heatmap_integ_pro.pdf"), height = 52, width = 80, units="mm", dpi = 300)

# Enhancers
triad_cors = triads %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("Enhancers")
ggsave(paste0(paths$out, "heatmap_integ_enh.pdf"), height = 52, width = 80, units="mm", dpi = 300)

##### DENSITIES: RNA VS ATAC #####

diads_ra = triads %>%
  dplyr::select(peakID:sig_rna) %>%
  distinct() %>%
  drop_na(corr_rna)
  
p1 = ggplot(diads_ra, aes(x=corr_atac, y=corr_rna))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("All")+
  theme(legend.position = "bottom")
p2 = diads_ra %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_atac, y=corr_rna))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Promoter")+
  theme(legend.position = "bottom")
p3 = diads_ra %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_atac, y=corr_rna))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Enhancer")+
  theme(legend.position = "bottom")
p1 + p2 + p3
ggsave(paste0(paths$out, "hist2d_rna_atac_age_corr.pdf"), height = 0.5*h, width = 1.5*w, units="mm", dpi = 300)

p1 = diads_ra %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_atac, fill=sig_rna))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_ra %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_atac, fill=sig_rna))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_atac_corr_in_rna_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

p1 = diads_ra %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_rna, fill=sig_atac))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_ra %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_rna, fill=sig_atac))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_rna_corr_in_atac_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

##### DENSITIES: ATAC VS MET #####

diads_am = triads %>%
  dplyr::select(probe_id, peakID, desc:element_type, corr_atac:sig_atac, corr_met:sig_met) %>%
  distinct() %>%
  drop_na(corr_met)
p1 = ggplot(diads_am, aes(x=corr_atac, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("All")+
  theme(legend.position = "bottom")
p2 = diads_am %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_atac, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Promoter")+
  theme(legend.position = "bottom")
p3 = diads_am %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_atac, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Enhancer")+
  theme(legend.position = "bottom")
p1 + p2 + p3
ggsave(paste0(paths$out, "hist2d_atac_met_age_corr.pdf"), height = 0.5*h, width = 1.5*w, units="mm", dpi = 300)

p1 = diads_am %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_atac, fill=sig_met))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_am %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_atac, fill=sig_met))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_atac_corr_in_met_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

p1 = diads_am %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_met, fill=sig_atac))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_am %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_met, fill=sig_atac))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_met_corr_in_atac_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

##### DENSITIES: RNA VS MET #####

diads_rm = triads %>%
  dplyr::select(probe_id, gene_ensemblID, element_type, corr_rna:sig_met) %>%
  distinct() %>%
  drop_na(corr_met, corr_rna)
p1 = ggplot(diads_rm, aes(x=corr_rna, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("All")+
  theme(legend.position = "bottom")
p2 = diads_rm %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_rna, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Promoter")+
  theme(legend.position = "bottom")
p3 = diads_rm %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_rna, y=corr_met))+
  geom_bin2d(bins=50)+
  geom_smooth(method = "lm", color="red")+
  stat_cor(color="red")+
  ggtitle("Enhancer")+
  theme(legend.position = "bottom")
p1 + p2 + p3
ggsave(paste0(paths$out, "hist2d_rna_met_age_corr.pdf"), height = 0.5*h, width = 1.5*w, units="mm", dpi = 300)

p1 = diads_rm %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_rna, fill=sig_met))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_rm %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_rna, fill=sig_met))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_rna_corr_in_met_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

p1 = diads_rm %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_met, fill=sig_rna))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Promoter")
p2 = diads_rm %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_met, fill=sig_rna))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c(cols$gray, cols$blue, cols$red))+
  ggtitle("Enhancer")
p1+p2+plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "hist_met_corr_in_rna_sig.pdf"), height = 51.8, width = 95.6, units="mm", dpi = 300)

##### KS TESTS: RNA VS ATAC #####

combos = expand.grid(c("rna", "atac"), c("rna", "atac"), c("Promoter", "Enhancer")) %>%
  dplyr::filter(Var1 != Var2) %>%
  arrange(Var1)
alt_map = list("Increasing" = "less", "Decreasing" = "greater")

res = data.frame()
for (i in 1:nrow(combos)) {
  v1 = combos[i, "Var1"]
  v2 = combos[i, "Var2"]
  e_type = combos[i, "Var3"]
  tmp = dplyr::filter(diads_ra, element_type == e_type)
  out = all_ks_tests(tmp, v1, v2, alt_map)
  out = mutate(out, element_type = e_type, .after = for_which_group)
  res = rbind(res, out)
}

cat("RNA VS ATAC KS TESTS\n")
res
cat("=====================\n")

##### KS TESTS: ATAC VS MET #####

combos = expand.grid(c("atac", "met"), c("atac", "met"), c("Promoter", "Enhancer")) %>%
  dplyr::filter(Var1 != Var2) %>%
  arrange(Var1)
alt_map = list("Increasing" = "greater", "Decreasing" = "less")

res = data.frame()
for (i in 1:nrow(combos)) {
  v1 = combos[i, "Var1"]
  v2 = combos[i, "Var2"]
  e_type = combos[i, "Var3"]
  tmp = dplyr::filter(diads_am, element_type == e_type)
  out = all_ks_tests(tmp, v1, v2, alt_map)
  out = mutate(out, element_type = e_type, .after = for_which_group)
  res = rbind(res, out)
}

cat("ATAC VS MET KS TESTS\n")
res
cat("=====================\n")

##### KS TESTS: RNA VS MET #####

combos = expand.grid(c("rna", "met"), c("rna", "met"), c("Promoter", "Enhancer")) %>%
  dplyr::filter(Var1 != Var2) %>%
  arrange(Var1)
alt_map = list("Increasing" = "greater", "Decreasing" = "less")

res = data.frame()
for (i in 1:nrow(combos)) {
  v1 = combos[i, "Var1"]
  v2 = combos[i, "Var2"]
  e_type = combos[i, "Var3"]
  tmp = dplyr::filter(diads_rm, element_type == e_type)
  out = all_ks_tests(tmp, v1, v2, alt_map)
  out = mutate(out, element_type = e_type, .after = for_which_group)
  res = rbind(res, out)
}

cat("RNA VS MET KS TESTS\n")
res
cat("=====================\n")

##### ATAC CLOCK VS HANNUM/HORVATH CLOCKS #####

# Hannum
hannum_overlap = triads %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(probe_id %in% rownames(hannum_coef))
hannum_overlap$coef_atac = final_coef[hannum_overlap$peakID, "X0"]
hannum_overlap$coef_hannum = hannum_coef[hannum_overlap$probe_id, "Coefficient"]
cat("Peaks shared with Hannum\n")
hannum_overlap
cat("=====================\n")

# Horvath
horvath_overlap = triads %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(probe_id %in% rownames(horvath_coef))
horvath_overlap$coef_atac = final_coef[horvath_overlap$peakID, "X0"]
horvath_overlap$coef_horvath = horvath_coef[horvath_overlap$probe_id, "CoefficientTraining"]
cat("Peaks shared with Horvath\n")
horvath_overlap
cat("=====================\n")

# Overlap of Hannum and Horvath
hh_overlap = triads %>%
  dplyr::filter(probe_id %in% rownames(hannum_coef)) %>%
  dplyr::filter(probe_id %in% rownames(horvath_coef))
hh_overlap$coef_hannum = hannum_coef[hh_overlap$probe_id, "Coefficient"]
hh_overlap$coef_horvath = horvath_coef[hh_overlap$probe_id, "CoefficientTraining"]
cat("Peaks shared between Hannum and Horvath\n")
hh_overlap
cat("=====================\n")

##### HEATMAPS: CLOCK SITES #####

triads_clock = triads %>%
  dplyr::filter(peakID %in% clock_all_peaks)

# Promoters
triad_cors = triads_clock %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("ATAC-clock promoters")
ggsave(paste0(paths$out, "heatmap_integ_pro_clock.pdf"), height = 52, width = 80, units="mm", dpi = 300)

# Enhancers
triad_cors = triads_clock %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("ATAC-clock enhancers")
ggsave(paste0(paths$out, "heatmap_integ_enh_clock.pdf"), height = 52, width = 80, units="mm", dpi = 300)

##### HEATMAPS: HANNUM CLOCK SITES #####

triads_hannum = triads %>%
  dplyr::filter(probe_id %in% rownames(hannum_coef))

# Promoters
triad_cors = triads_hannum %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("Hannum promoters")
ggsave(paste0(paths$out, "heatmap_integ_pro_hannum.pdf"), height = 52, width = 80, units="mm", dpi = 300)

# Enhancers
triad_cors = triads_hannum %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
triad_cors$cors %>% 
  rownames_to_column(var="Feature1") %>%
  pivot_longer(cols=-Feature1, names_to = "Feature2") %>%
  mutate(Feature1 = factor(Feature1, levels=c("r(RNA-age)", "r(ATAC-age)", "r(MET-age)"))) %>%
  mutate(Feature2 = factor(Feature2, levels=c("r(MET-age)", "r(ATAC-age)", "r(RNA-age)"))) %>%
  ggplot(., aes(x=Feature1, y=Feature2, fill=value, label = round(value, 3)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low=cols$blue, mid=cols$gray, high=cols$red, limits=c(-1,1))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("Hannum enhancers")
ggsave(paste0(paths$out, "heatmap_integ_enh_hannum.pdf"), height = 52, width = 80, units="mm", dpi = 300)

##### SCATTERS #####

p1 = diads_ra %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_atac, y=corr_rna))+
  geom_point(aes(color=sig_atac))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p2 = diads_rm %>%
  merge(., cpg_peak_pairs[c("probe_id", "peak_id")], by="probe_id") %>%
  dplyr::filter(peak_id %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_met, y=corr_rna))+
  geom_point(aes(color=sig_met))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p3 = diads_am %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Promoter") %>%
  ggplot(., aes(x=corr_met, y=corr_atac))+
  geom_point(aes(color=sig_met))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p1+p2+p3+ plot_annotation(title = "Clock Promoters")
ggsave(paste0(paths$out, "corr_clock_integ_pro.pdf"), height = 60.9, width = w, units="mm", dpi = 300)


p1 = diads_ra %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_atac, y=corr_rna))+
  geom_point(aes(color=sig_atac))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p2 = diads_rm %>%
  merge(., cpg_peak_pairs[c("probe_id", "peak_id")], by="probe_id") %>%
  dplyr::filter(peak_id %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_met, y=corr_rna))+
  geom_point(aes(color=sig_met))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p3 = diads_am %>%
  dplyr::filter(peakID %in% clock_all_peaks) %>%
  dplyr::filter(element_type == "Enhancer") %>%
  ggplot(., aes(x=corr_met, y=corr_atac))+
  geom_point(aes(color=sig_met))+
  scale_color_manual(values=c(cols$dark_gray, cols$blue, cols$red))+
  geom_smooth(method="lm", color="black")+
  stat_cor(label.y = 0.7, cor.coef.name="r")+
  guides(color="none")
p1+p2+p3+ plot_annotation(title = "Clock Enhancers")
ggsave(paste0(paths$out, "corr_clock_integ_enh.pdf"), height = 60.9, width = w, units="mm", dpi = 300)

closeAllConnections()
