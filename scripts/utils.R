# Convenience function to make contingency table
make_contingency = function(hit_in_fg, fg, hit_in_bg, bg) {
  n_fg = length(fg)
  n_hit_fg = length(hit_in_fg)
  n_miss_fg = n_fg - n_hit_fg
  
  rest = setdiff(bg, fg)
  hit_in_rest = setdiff(hit_in_bg, hit_in_fg)
  n_rest = length(rest)
  n_hit_rest = length(hit_in_rest)
  n_miss_rest = n_rest - n_hit_rest

  return(matrix(c(n_hit_fg, n_miss_fg, n_hit_rest, n_miss_rest), nrow=2))
}

# Group peaks into categories (promoter, enhancer ...)
# Creates lists of each category, also counts them for convenience
element_type_lists = function(peaks, anno) {
  anno = anno[!is.na(anno$element_type), ]
  anno = anno[anno$peakID %in% peaks,]
  # Peak sublists for various element types
  promoter_peaks = unique(anno$peakID[anno$element_type == "Promoter"])
  enhancer_peaks = unique(anno$peakID[anno$element_type == "Enhancer"])
  pro_and_enh_peaks = intersect(promoter_peaks, enhancer_peaks)
  promoter_only_peaks = setdiff(promoter_peaks, pro_and_enh_peaks)
  enhancer_only_peaks = setdiff(enhancer_peaks, pro_and_enh_peaks)
  unanno_peaks = setdiff(peaks, union(promoter_peaks, enhancer_peaks))
  res = list("all_peaks" = peaks,
             "all_peaks_n" = length(peaks),
             "pro_peaks" = promoter_peaks,
             "pro_peaks_n" = length(promoter_peaks),
             "enh_peaks" = enhancer_peaks,
             "enh_peaks_n" = length(enhancer_peaks),
             "pro_and_enh_peaks" = pro_and_enh_peaks,
             "pro_and_enh_peaks_n" = length(pro_and_enh_peaks),
             "pro_only_peaks" = promoter_only_peaks,
             "pro_only_peaks_n" = length(promoter_only_peaks),
             "enh_only_peaks" = enhancer_only_peaks,
             "enh_only_peaks_n" = length(enhancer_only_peaks),
             "unanno_peaks" = unanno_peaks,
             "unanno_peaks_n" = length(unanno_peaks))
        
  return(res)
}

# Get all genes linked to list of peaks, provided annotation table
peaks_to_genes = function(peaks, anno, ensembl = FALSE) {
  if (ensembl) {
    genes = unique(anno$gene_ensemblID[anno$peakID %in% peaks])
  } else {
    genes = unique(anno$gene_symbol[anno$peakID %in% peaks])
  }
  genes = genes[!is.na(genes)]
  return(genes)
}

# Get all peaks linked to list of genes, provided annotation table
genes_to_peaks = function(genes, anno) {
  peaks = unique(anno$peakID[anno$gene_symbol %in% genes])
  peaks = peaks[!is.na(peaks)]
  return(peaks)
}

# Correlation tests of all matrix cols against one variable
cor_test_matrix = function(matrix, target, type) {
  r = c()
  pvals = c()
  for (i in 1:ncol(matrix)) {
    var = matrix[,i]
    test = cor.test(var, target, method=type)
    r = append(r, test$estimate)
    pvals = append(pvals, test$p.value)
  }
  return(data.frame(r, pvals, row.names = colnames(matrix)))
}

# Creates list of overlapping peak-probe pairs
# Requires dfs of peaks with id, start, end columns
#              of probes with id, start
# This assumes coordinared to refer to the same chromosome
map_cpg_to_peak = function(peaks, probes) {
  mappings = list()
  j=1
  for (i in 1:nrow(peaks)) {
    peakId = rownames(peaks)[i]
    peak_start = peaks$start[i]
    peak_end = peaks$end[i]
    while (j < nrow(probes)) {
      probeId = rownames(probes)[j]
      probe_pos = probes$start[j]
      if (probe_pos > peak_end) {
        break
      } else if (probe_pos >= peak_start) {
        mappings = append(mappings, list(c(peakId, probeId)))
      }
      j = j+1
    }
  }
  if (length(mappings) != 0) {
    mappings = data.frame(t(do.call(cbind.data.frame, mappings)))
    colnames(mappings) = c("PeakId", "ProbeId")
    rownames(mappings) = 1:nrow(mappings)
    return(mappings)
  } else {
    return(NULL)
  }
}

triad_heatmap = function(triads) {
  cors = data.frame()
  pvals = data.frame()
  
  # Diagonal
  cors["r(RNA-age)", "r(RNA-age)"] = 1
  cors["r(ATAC-age)", "r(ATAC-age)"] = 1
  cors["r(MET-age)", "r(MET-age)"] = 1
  pvals["r(RNA-age)", "r(RNA-age)"] = 0
  pvals["r(ATAC-age)", "r(ATAC-age)"] = 0
  pvals["r(MET-age)", "r(MET-age)"] = 0
  
  # RNA VS ATAC
  tmp = triads %>%
    dplyr::select(gene_ensemblID, peakID, corr_rna, corr_atac) %>%
    dplyr::distinct() %>%
    drop_na()
  test = cor.test(tmp$corr_rna, tmp$corr_atac)
  cors["r(RNA-age)", "r(ATAC-age)"] = test$estimate
  cors["r(ATAC-age)", "r(RNA-age)"] = test$estimate
  pvals["r(RNA-age)", "r(ATAC-age)"] = test$p.value
  pvals["r(ATAC-age)", "r(RNA-age)"] = test$p.value
  
  # MET VS ATAC
  tmp = triads %>%
    dplyr::select(probe_id, peakID, corr_met, corr_atac) %>%
    dplyr::distinct() %>%
    drop_na()
  test = cor.test(tmp$corr_met, tmp$corr_atac)
  cors["r(MET-age)", "r(ATAC-age)"] = test$estimate
  cors["r(ATAC-age)", "r(MET-age)"] = test$estimate
  pvals["r(MET-age)", "r(ATAC-age)"] = test$p.value
  pvals["r(ATAC-age)", "r(MET-age)"] = test$p.value
  
  # MET VS RNA
  tmp = triads %>%
    dplyr::select(probe_id, gene_ensemblID, corr_met, corr_rna) %>%
    dplyr::distinct() %>%
    drop_na()
  test = cor.test(tmp$corr_met, tmp$corr_rna)
  cors["r(MET-age)", "r(RNA-age)"] = test$estimate
  cors["r(RNA-age)", "r(MET-age)"] = test$estimate
  pvals["r(MET-age)", "r(RNA-age)"] = test$p.value
  pvals["r(RNA-age)", "r(MET-age)"] = test$p.value
  
  res = list()
  res$cors = cors
  res$pvals = pvals
  return(res)
}

all_ks_tests = function(diads, v1, v2, alt_map) {
  KS_res = data.frame()
  corr = paste0("corr_", v1)
  sig = paste0("sig_", v2)
  # INCREASING
  test = ks.test(dplyr::filter(diads, get(sig) == "Increasing")[,corr],
                 dplyr::filter(diads, get(sig) == "Constant")[,corr],
                 alternative = alt_map[["Increasing"]])
  KS_res[1, "which_cor"] = corr
  KS_res[1, "for_which_group"] = paste0(sig, "_increasing")
  KS_res[1, "alt"] = alt_map[["Increasing"]]
  KS_res[1, "D"] = test$statistic
  KS_res[1, "pval"] = test$p.value
  # DECREASING
  test = ks.test(dplyr::filter(diads, get(sig) == "Decreasing")[,corr],
                 dplyr::filter(diads, get(sig) == "Constant")[,corr],
                 alternative = alt_map[["Decreasing"]])
  KS_res[2, "which_cor"] = corr
  KS_res[2, "for_which_group"] = paste0(sig, "_decreasing")
  KS_res[2, "alt"] = alt_map[["Decreasing"]]
  KS_res[2, "D"] = test$statistic
  KS_res[2, "pval"] = test$p.value
  return(KS_res)
}