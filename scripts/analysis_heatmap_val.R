library(matrixStats)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(edgeR)

setwd("/scratch/fmorandi/ChromAcc-clock")
source("./scripts/utils.R")

paths = list()
paths$data_us = "./data/paper_data/"
paths$rna_addon = "./data/paper_addon_rna/"
paths$rna_ucar = "./data/pipeline_outputs/rna_ucar/03_counts/counts.tsv"
paths$meta_ucar = "./data/meta_ucar.tsv"
paths$out = "./data/paper_data/outputs_2023-03-20_13-43-53/"

enhancer_links = "closest"

# RNA data from 10k immunomes, RNA-seq
rna_alt1 = fread(paste0(paths$rna_addon, "pbmc_rna_counts.csv"))
rna_alt1 = data.frame(rna_alt1, row.names = 1)
dge = DGEList(counts=rna_alt1)
dge = calcNormFactors(dge)
to_keep = filterByExpr(dge)
dge = dge[to_keep, ]
rna_alt1 = cpm(dge, normalized.lib.sizes=T)
rna_alt1 = data.frame(t(rna_alt1)) # <!> not log normalized but doesnt matter for Spearman
rm(dge, to_keep)

# RNA data from 10k immunomes, Affymetrix
rna_alt2 = fread(paste0(paths$rna_addon, "pbmc_affy_normalized.csv"))
rna_alt2 = data.frame(rna_alt2, row.names = 2)
meta_alt2 = rna_alt2[, 1:7]
rna_alt2 = rna_alt2[, -c(1:7)]

# RNA data from ucar
rna_alt3 = fread(paths$rna_ucar)
rna_alt3 = data.frame(rna_alt3, row.names = 1)
rna_alt3 = rna_alt3[-c(1:5)]
colnames(rna_alt3) = str_extract(colnames(rna_alt3), "(SRR[0-9]+)", group=1)
dge = DGEList(counts=rna_alt3)
dge = calcNormFactors(dge)
to_keep = filterByExpr(dge)
dge = dge[to_keep, ]
rna_alt3 = cpm(dge, normalized.lib.sizes=T)
rna_alt3 = data.frame(t(rna_alt3)) # <!> not log normalized but doesnt matter for Spearman
rm(dge, to_keep)

# Our precomputed ATAC Spearman cors
peak_info = read.table(paste0(paths$data_us, "outputs_2023-03-11_13-34-19/peak_info_cors.tsv"), header=T)

# Alt rna meta1
meta_alt1 = read.csv(paste0(paths$rna_addon, "pbmc_rna_meta.csv"))
# meta_alt2 = read.csv(paste0(paths$rna_addon, "pbmc_affy_meta.csv")) # Let's use the one included with the array
meta_alt3 = read.table(paths$meta_ucar, sep="\t", header = T, row.names = 1) %>%
  dplyr::filter(ANALYTE_TYPE == "RNA-seq") %>%
  mutate(Age = as.numeric(str_replace_all(Age, "\\+", "")))

##### ALT RNA SPEARMAN CORRELATION AND SAVE #####

sig_th_rna = 0.01

all(rownames(rna_alt1) == meta_alt1$Subject)
# all(rownames(rna_alt2) == meta_alt1$Subject) # Guaranteed because they were on the same table
rna_alt3 = rna_alt3[meta_alt3$SRR, ]
all(rownames(rna_alt3) == meta_alt3$SRR)

rna_age_cors1 = cor_test_matrix(rna_alt1, meta_alt1$Age, "spearman")
rna_age_cors2 = cor_test_matrix(rna_alt2, meta_alt2$age, "spearman")
rna_age_cors3 = cor_test_matrix(rna_alt3, meta_alt3$Age, "spearman")

rna_age_cors1$padj = p.adjust(rna_age_cors1$pvals, method = "fdr")
rna_age_cors2$padj = p.adjust(rna_age_cors2$pvals, method = "fdr")
rna_age_cors3$padj = p.adjust(rna_age_cors3$pvals, method = "fdr")

colnames(rna_age_cors1) = c("corr","pval","padj")
colnames(rna_age_cors2) = c("corr","pval","padj")
colnames(rna_age_cors3) = c("corr","pval","padj")

write.table(rna_age_cors1, paste0(paths$rna_addon, "rna_age_cors_1.tsv"))
write.table(rna_age_cors2, paste0(paths$rna_addon, "rna_age_cors_2.tsv"))
write.table(rna_age_cors3, paste0(paths$rna_addon, "rna_age_cors_3.tsv"))

##### LOAD ADDITIONAL FILES NEEDED FOR HEATMAP #####

rna_age_cors0 = read.table(paste0(paths$data_us, "outputs_2023-03-16_10-45-59/gene_info_cors.tsv"), row.names = 2)
rna_age_cors1 = read.table(paste0(paths$rna_addon, "rna_age_cors_1.tsv"))
rna_age_cors2 = read.table(paste0(paths$rna_addon, "rna_age_cors_2.tsv"))
rna_age_cors3 = read.table(paste0(paths$rna_addon, "rna_age_cors_3.tsv"))

hannum_age_cors = read.table(paste0(paths$data_us, "hannum_age_cors.tsv"), header = 1)
cpg_peak_pairs = read.table(paste0(paths$data_us, "peak_cpg_map.tsv"), header = 1)
anno = read.table(paste0(paths$data_us, "peak_anno_", enhancer_links, ".txt"), header=T)

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

##### PREP TRIADS #####

length(intersect(rownames(rna_age_cors0), anno$gene_ensemblID))
length(intersect(rownames(rna_age_cors1), anno$gene_symbol))
length(intersect(rownames(rna_age_cors2), anno$gene_symbol))
length(intersect(rownames(rna_age_cors3), anno$gene_ensemblID))

prep_triads = function(anno, peak_info, gene_info, met_info, met_links, gene_ids="gene_symbol") {
  triads = merge(anno, peak_info[,c("corr", "pval", "padj")], by.x = "peakID", by.y = "row.names")
  triads = merge(triads, gene_info, by.x = gene_ids, by.y = 0,
                 suffixes = c("_atac", "_rna"), all.x=T)
  triads = merge(triads, met_links, by.x = "peakID", by.y = "peak_id", all.x=T)
  triads = met_info %>%
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
  return(triads)
}

triads0 = prep_triads(anno, peak_info, rna_age_cors0, hannum_age_cors, cpg_peak_pairs, gene_ids = "gene_ensemblID")
triads1 = prep_triads(anno, peak_info, rna_age_cors1, hannum_age_cors, cpg_peak_pairs, gene_ids = "gene_symbol")
triads2 = prep_triads(anno, peak_info, rna_age_cors2, hannum_age_cors, cpg_peak_pairs, gene_ids = "gene_symbol")
triads3 = prep_triads(anno, peak_info, rna_age_cors3, hannum_age_cors, cpg_peak_pairs, gene_ids = "gene_ensemblID")

##### HEATMAP 0 #####

triad_cors = triads0 %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
hp0 = triad_cors$cors %>% 
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

# Enhancers
triad_cors = triads0 %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
he0 = triad_cors$cors %>% 
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

##### HEATMAP 1 #####

triad_cors = triads1 %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
hp1 = triad_cors$cors %>% 
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

# Enhancers
triad_cors = triads1 %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
he1 = triad_cors$cors %>% 
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

##### HEATMAP 2 #####

triad_cors = triads2 %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
hp2 = triad_cors$cors %>% 
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

# Enhancers
triad_cors = triads2 %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
he2 = triad_cors$cors %>% 
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

##### HEATMAP 3 #####

triad_cors = triads3 %>%
  dplyr::filter(element_type == "Promoter") %>%
  triad_heatmap(.)
hp3 = triad_cors$cors %>% 
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

# Enhancers
triad_cors = triads3 %>%
  dplyr::filter(element_type == "Enhancer") %>%
  triad_heatmap(.)
he3 = triad_cors$cors %>% 
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

##### FINAL PLOT #####

hp0 + hp1 + hp2 + hp3 + plot_layout(nrow=2, guides = "collect")
he0 + he1 + he2 + he3 + plot_layout(nrow=2, guides = "collect")

# tl = hp3 + theme(axis.text.x = element_blank())
# tr = he3 + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
# bl = hp1 + theme(plot.title = element_blank())
# br = he1 + theme(plot.title = element_blank(), axis.text.y = element_blank())
# tl+tr+bl+br+plot_layout(nrow=2, guides = "collect")
top = hp3 + theme(axis.text.x = element_blank())
bot = he3
top / bot  + plot_layout(guides="collect", ) & theme(legend.position = "bottom")
ggsave(paste0(paths$out, "heatmap_integ_alt_rna.pdf"), height = 100, width = 60, units="mm", dpi = 300)
