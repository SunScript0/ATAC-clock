library(tidyverse)
library(edgeR)
library(sp)
library(patchwork)

setwd("/scratch/fmorandi/ChromAcc-clock")

paths = list()
paths$pipeline ="./data/pipeline_outputs/rna/"
paths$out="./data/paper_data/"

# Options
w = 174 # mm
h = 230

dir.create(paste0(paths$out, "figures"), showWarnings = F)

##### LOAD DATA #####

counts = read.csv(paste0(paths$pipeline, "GSE193141_rna_counts_011622.csv"), header=T)
gene_info = counts[, 1:2]
counts = counts[, -c(1:2)]
colnames(counts) = str_replace_all(colnames(counts), "RNA_", "")
rownames(gene_info) = gene_info$Geneid
rownames(counts) = rownames(gene_info)
counts = data.frame(t(counts))
# Already remove genes that are zero throughout
counts = counts[, colSums(counts != 0) > 0]
gene_info = gene_info[colnames(counts), ]

meta = read.table(paste0(paths$out, "meta_atac.tsv"), header=T)

##### NORMALIZE #####

dge = DGEList(counts=t(counts))
dge = calcNormFactors(dge)

tmm = cpm(dge, normalized.lib.sizes=T)
tmm = data.frame(t(tmm))

##### FILTER GENES #####

to_keep = filterByExpr(dge)
tmm = tmm[, to_keep]
gene_info = gene_info[to_keep, ]

##### REMOVE OUTLIERS #####

pca = prcomp(log1p(tmm))
pca = merge(pca$x, meta, by.x=0, by.y="Subject")
p = ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point()+
  stat_ellipse(level=0.99)
ellipse = ggplot_build(p)$data[[2]]
pca$PassesQC_rna = as.logical(point.in.polygon(pca$PC1, pca$PC2, ellipse$x, ellipse$y))
ggplot(pca, aes(x=PC1, y=PC2))+
  geom_point(aes(color=PassesQC_rna))+
  stat_ellipse(level=0.99)+
  theme_light()
ggsave(paste0(paths$out, "figures/qc_rna_ellipse.pdf"), height = 0.4*h, width = 1*w, units="mm", dpi = 300)

##### SUMMARIZE CLEANUP #####

all(meta$Subject == pca$Row.names)
meta$PassesQC_rna = pca$PassesQC_rna

p1 = ggplot(meta, aes(x=Age, fill=!PassesQC_rna))+
  geom_histogram(breaks=0.5+seq(20, 74))+
  labs(fill="Failing QC") +
  theme_light()
p2 = ggplot(meta, aes(x=Sex, fill=!PassesQC_rna))+
  geom_bar(position="fill")+
  labs(fill="Failing QC") +
  theme_light()
(p1 | p2) + plot_layout(widths = c(4, 1), guides = 'collect')
ggsave(paste0(paths$out, "figures/qc_rna_distrs.pdf"), height = 0.3*h, width = w, units="mm", dpi = 300)

meta$PassesQC = meta$PassesQC_atac & meta$PassesQC_rna

p1 = ggplot(meta, aes(x=Age, fill=!PassesQC))+
  geom_histogram(breaks=0.5+seq(20, 74))+
  labs(fill="Failing QC") +
  theme_light()
p2 = ggplot(meta, aes(x=Sex, fill=!PassesQC))+
  geom_bar(position="fill")+
  labs(fill="Failing QC") +
  theme_light()
(p1 | p2) + plot_layout(widths = c(4, 1), guides = 'collect')
ggsave(paste0(paths$out, "figures/qc_both_distrs.pdf"), height = 0.3*h, width = w, units="mm", dpi = 300)

##### SAVE OUTPUTS #####

writeLines(meta[!meta$PassesQC_rna, "Subject"], paste0(paths$out, "qc_failed_rna.txt"))

write.table(tmm, file=paste0(paths$out, "data_rna_tmm.tsv"), quote=F, sep="\t")
write.table(gene_info, file=paste0(paths$out, "gene_info.tsv"), quote=F, sep="\t")
write.table(meta, file=paste0(paths$out, "meta_final.tsv"), quote=F, sep="\t")
