# Helper to build ranked list and run fgsea
run_gsea <- function(tt, pathways, gene_col = "SYMBOL", nPermSimple = 10000) {
  tt <- tt[!is.na(tt[[gene_col]]) & tt[[gene_col]] != "", ]
  
  n_before <- nrow(tt)
  tt <- tt[order(-tt$AveExpr), ]
  tt <- tt[!duplicated(tt[[gene_col]]), ]
  message("Removed ", n_before - nrow(tt), " duplicate gene names (kept highest AveExpr)")
  
  ranks <- setNames(tt$t, tt[[gene_col]])
  ranks <- sort(ranks, decreasing = TRUE)
  
  res <- fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
  res <- res[order(res$pval), ]
  return(res)
}