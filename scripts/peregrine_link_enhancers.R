library(plyr)
library(tidyverse)
library(gtools)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
# library(biomaRt)

setwd("/scratch/fmorandi/ChromAcc-clock")

##### SETTINGS #####

blood_links_only = F
perepath = "./data/external_resources/peregrine/" # Peregrine data path
datapath = "./data/paper_data/" # Pipeline output folder

##### LOAD DATA #####

# Peregrine
enhancers = read.table(paste0(perepath, "PEREGRINEenhancershg38"), header=F)
colnames(enhancers) = c("chr", "start", "end", "enhancerID")
links = read.table(paste0(perepath, "enh_gene_link_assay_tissue_pval_snp_score"), header=T, fill=T, sep="\t")
tissues = read.table(paste0(perepath, "PEREGRINEtissues"), header=F, sep = "\t")
if (blood_links_only) {
  # Tissue 48 is blood
  links = links[links$tissue == 48, ]
}

# Peak data
peaks = read.table(paste0(datapath, "peak_info.tsv"))
colnames(peaks) = c("desc", "chr", "start", "end", "strand", "len")
peaks$peakID = rownames(peaks)

##### LINK PEAKS, ENHANCERS AND GENES #####

peak_ranges = GRanges(
  seqnames = peaks$chr,
  ranges = IRanges(peaks$start, end = peaks$end, names = peaks$peakID))
enh_ranges = GRanges(
  seqnames = enhancers$chr,
  ranges = IRanges(enhancers$start, end = enhancers$end, names = enhancers$enhancerID))

mappings = data.frame(findOverlaps(peak_ranges, enh_ranges))
colnames(mappings) = c("peakID", "enhancerID")
mappings$peakID=names(peak_ranges)[mappings$peakID]
mappings$enhancerID=names(enh_ranges)[mappings$enhancerID]

compact_links = links %>%
  dplyr::select(enhancer, gene, tissue) %>%
  dplyr::group_by(enhancer, gene) %>%
  dplyr::summarize(val_in_blood = 48 %in% unique(tissue))
    
mappings = merge(mappings, compact_links, by.x = "enhancerID", 
                 by.y = "enhancer", all.x = T, all.y = F, sort = F)
mappings = merge(peaks, mappings, by="peakID", all=T)

rm(compact_links, enh_ranges, enhancers, links, tissues)

##### CHIP SEEKER #####

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
anno = annotatePeak(peak_ranges, 
                    TxDb = txdb,
                    annoDb = "org.Hs.eg.db",
                    tssRegion = c(-1000, 1000))
plotAnnoPie(anno)
annoDF = as.data.frame(anno)

# merge annotations to list of peaks by region to recover peakID
annoDF = merge(peaks, annoDF, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"))
annoDF = annoDF[c("peakID", "chr", "start", "end", "len", "desc", "annotation", "geneId", "ENSEMBL", "SYMBOL")]
colnames(annoDF) = c("peakID", "chr", "start", "end", "len", "desc", "element_type", "gene_entrezID", "gene_ensemblID", "gene_symbol")

##### MERGE CHIP SEEKER AND ENHANCER LINKS #####

# Some peaks have enhancer annotation but no gene link:
sum(is.na(mappings$enhancerID))
sum(is.na(mappings$gene))

# Restrict to peaks which have enhancer annotations
enhancer_links = mappings[!is.na(mappings$enhancerID), c("peakID", "chr", "start", "end", "enhancerID", "val_in_blood")]
# Split to get gene ids
long_geneIDs = data.frame(str_split_fixed(mappings$gene, "\\|", 3))
long_geneIDs = long_geneIDs[!is.na(mappings$enhancerID), ]
colnames(long_geneIDs) = c("Species", "gene_hgncID", "gene_uniprotID")
enhancer_links = cbind(enhancer_links, long_geneIDs[c("gene_hgncID", "gene_uniprotID")])
enhancer_links[enhancer_links == ""] = NA
rm("long_geneIDs")

promoter_peaks = annoDF$peakID[grep("promoter", annoDF$element_type, ignore.case = TRUE)]
enhancer_peaks = unique(enhancer_links$peakID)
pro_enh_peaks = intersect(promoter_peaks, enhancer_peaks)
unannotated_peaks = setdiff(peaks$peakID, union(promoter_peaks, enhancer_peaks))

pro_enh_n = length(pro_enh_peaks)
pro_n = length(promoter_peaks) - pro_enh_n
enh_n = length(enhancer_peaks) - pro_enh_n
unanno_n = length(unannotated_peaks)

pie(c(pro_n, pro_enh_n, enh_n, unanno_n),
    c("Promoter", "Promoter and enhancer", "Enhancer", "Unannotated"))

annoDF_promoter_only = annoDF[grep("promoter", annoDF$element_type, ignore.case = TRUE),]

enhancer_links$element_type = "Enhancer"

merged_table = rbind.fill(enhancer_links, annoDF_promoter_only) 
# Here some enhancers have a row with many NAs (unlinked enhancer) but also links
# Theses arent really needed since other rows already say they are enhancers
nrow(merged_table) == nrow(enhancer_links) + nrow(annoDF_promoter_only)
# merged_table = merged_table[mixedorder(merged_table$peakID), ]
merged_table = merged_table %>%
  dplyr::select(-enhancerID) %>% 
  distinct() # drop enhancer ID, makes some rows redundant

##### CONVERT GENE IDS #####

merged_table$gene_hgncID = str_replace_all(merged_table$gene_hgncID, "HGNC=", "")
merged_table$gene_uniprotID = str_replace_all(merged_table$gene_uniprotID, "UniProtKB=", "")

# Note: conversions dont always agree!
symbols = mapIds(org.Hs.eg.db, keys = merged_table$gene_uniprotID, keytype = "UNIPROT", column = 'ENSEMBL')
merged_table[!is.na(merged_table$gene_uniprotID), "gene_ensemblID"] = as.character(symbols[na.omit(merged_table$gene_uniprotID)])

symbols = mapIds(org.Hs.eg.db, keys = merged_table$gene_uniprotID, keytype = "UNIPROT", column = 'SYMBOL')
merged_table[!is.na(merged_table$gene_uniprotID), "gene_symbol"] = as.character(symbols[na.omit(merged_table$gene_uniprotID)])

symbols = mapIds(org.Hs.eg.db, keys = merged_table$gene_uniprotID, keytype = "UNIPROT", column = 'ENTREZID')
merged_table[!is.na(merged_table$gene_uniprotID), "gene_entrezID"] = as.character(symbols[na.omit(merged_table$gene_uniprotID)])

# mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# symbols2 = getBM(
#   attributes = c("ensembl_gene_id","uniprot_gn_id"), 
#   filters = "uniprot_gn_id", 
#   values = merged_table$gene_uniprotID, 
#   mart = mart)

##### CLEANUP AND SAVE #####

# Unannotated peaks are missing from merged table
nrow(peaks) == length(unique(merged_table$peakID)) + unanno_n

merged_table = merge(peaks, merged_table[, c("peakID", "element_type", "val_in_blood", 
                                              "gene_entrezID", "gene_ensemblID", "gene_symbol")], 
                      all.x=TRUE) %>%
  dplyr::arrange(peakID, element_type, val_in_blood) %>%
  mutate(dupe = duplicated(paste(peakID, element_type))) %>%
  dplyr::filter(!(dupe & is.na(val_in_blood))) %>% # Remove rows with NA link if they already have other links
  dplyr::select(-dupe)

merged_table = merged_table[mixedorder(merged_table$peakID), ]

length(unique(merged_table$peakID))
sum(is.na(merged_table$element_type))
nrow(distinct(merged_table))

write.table(merged_table, file=paste0(datapath, "/peak_anno_peregrine.txt"), row.names = FALSE)

# If !val_in_blood remove gene links and save a copy
merged_table_blood  = merged_table
merged_table_blood$gene_entrezID[merged_table_blood$val_in_blood==F] = NA
merged_table_blood$gene_ensemblID[merged_table_blood$val_in_blood==F] = NA
merged_table_blood$gene_symbol[merged_table_blood$val_in_blood==F] = NA
merged_table_blood = distinct(merged_table_blood)
length(unique(merged_table_blood$peakID))

write.table(merged_table_blood, file=paste0(datapath, "/peak_anno_blood.txt"), row.names = FALSE)

# Simply link to closest gene
merged_table2 = merged_table_blood %>%
  dplyr::select(peakID:element_type) %>%
  distinct() %>%
  left_join(annoDF[, c("peakID", "gene_entrezID", "gene_ensemblID", "gene_symbol")], by="peakID")
length(unique(merged_table2$peakID))

write.table(merged_table2, file=paste0(datapath, "/peak_anno_closest.txt"), row.names = FALSE)
