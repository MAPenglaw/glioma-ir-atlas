#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(data.table)
})

## ============================================================
## G4：GBM 内部独立做 GO / KEGG 富集（摆脱乳腺癌 GSE227187 的唯一注释来源）
##
## 三个排序式 GSEA（基因 ranking → GSEA），全部只用 TCGA-GBM 数据：
##   A. retention      : 基因的 IR 事件 meanPSI 中位数（滞留倾向）
##   B. variability    : 基因的 IR 事件 stdPSI   中位数（跨肿瘤可塑性）
##   C. nmd_suscept    : 基因的 PSI–mRNA 耦合 rho（越负越靠前，主）
##                       + 敏感性分析用 rho_partial（已扣除文库复杂度）
## ============================================================

base <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
OUT  <- file.path(base, "ri_atlas/output")
FIG  <- file.path(base, "ri_atlas/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SRC <- paste0("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/",
              "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")

cat("[1] 读取 PSI meta ...\n")
psi <- suppressWarnings(fread(cmd = paste0("zcat < '", SRC, "/SplAdder_TCGA_GBM_PSI.txt.gz'"),
                              sep = "\t", header = TRUE, nThread = 2,
                              select = c("GeneSymbol", "SpliceType", "stdPSI", "MedianPSI", "pctSamples")))
ir  <- psi[psi$SpliceType == "IR", ]
ir  <- ir[!is.na(MedianPSI) & !is.na(stdPSI) & !is.na(GeneSymbol) & GeneSymbol != "", ]
cat("    IR 事件：", nrow(ir), "｜宿主基因：", uniqueN(ir$GeneSymbol), "\n")

coup <- fread(file.path(OUT, "G1_psi_expr_coupling.tsv"))
cat("    G1 耦合事件：", nrow(coup), "\n")

## ---------- 2. 基因级排序 --------------------------------------------
cat("[2] 基因级聚合 → 三条 ranking\n")
agg <- ir[, .(nIR = .N,
              psi_med   = median(MedianPSI),
              psi_max   = max(MedianPSI),
              sd_med    = median(stdPSI)), by = GeneSymbol]

coup_g <- coup[!is.na(gene) & gene != "", .(gene = gene,
                                            rho = rho_psi_expr,
                                            rho_partial = rho_partial)]
## 一个基因可能多个事件 → 取中位数 rho
coup_agg <- coup_g[, .(rho = median(rho, na.rm = TRUE),
                       rho_partial = median(rho_partial, na.rm = TRUE)), by = gene]

setkey(agg, GeneSymbol); setkey(coup_agg, gene)
DT <- merge(agg, coup_agg, by.x = "GeneSymbol", by.y = "gene", all.x = TRUE)

## 只在至少 1 个 IR 事件且表达可用的基因上做（universe）
cat("    合并后基因：", nrow(DT), "｜有耦合值：", sum(!is.na(DT$rho)), "\n")

mk_rank <- function(x) {
  x <- sort(x, decreasing = TRUE)
  return(x)
}
rankings <- list(
  retention   = setNames(DT$psi_med,  DT$GeneSymbol),
  variability = setNames(DT$sd_med,   DT$GeneSymbol),
  nmd_suscept = setNames(-DT$rho,     DT$GeneSymbol),   # 越负越靠前
  nmd_partial = setNames(-DT$rho_partial, DT$GeneSymbol)
)
rankings <- lapply(rankings, function(v) { v <- v[!is.na(v)]; sort(v, decreasing = TRUE) })
cat("    ranking 长度：", sapply(rankings, length), "\n")

## ---------- 3. GSEA -------------------------------------------------

run_gsea <- function(geneList, ont = "BP", pvalueCutoff = 0.05, Max = 500) {
  gseGO(geneList       = geneList,
        ont            = ont,
        OrgDb          = org.Hs.eg.db,
        keyType        = "SYMBOL",
        pvalueCutoff   = pvalueCutoff,
        pAdjustMethod  = "BH",
        minGSSize      = 20,
        maxGSSize      = Max,
        eps            = 0,
        verbose        = FALSE)
}

res <- list()
cat("\n[3] 运行 GSEA (GO-BP / GO-CC / GO-MF)\n")
for (nm in names(rankings)) {
  for (ont in c("BP", "CC", "MF")) {
    cat("    ", nm, "-", ont, " ... ", sep = "")
    r <- try(run_gsea(rankings[[nm]], ont), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r) || nrow(r) == 0) { cat("无结果\n"); next }
    cat(nrow(r), "条\n")
    res[[paste(nm, ont, sep = "_")]] <- r
    r@result$signature <- nm
    write.table(r@result, file.path(OUT, sprintf("G4_GSEA_%s_%s.tsv", nm, ont)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
cat("    完成 ", length(res), " 组\n", sep = "")

## ---------- 4. 简化展示 ---------------------------------------------
cat("\n[4] 各 signature 的 top 通路（按 |NES| 排序，取前 8）\n")
topall <- list()
for (k in names(res)) {
  r <- res[[k]]@result
  r <- r[order(abs(r$NES), decreasing = TRUE), ][1:min(8, nrow(r)), ]
  cat("\n=== ", k, " ===\n", sep = "")
  print(r[, c("Description", "setSize", "NES", "pvalue", "p.adjust")], row.names = FALSE)
  r$signature <- k
  topall[[k]] <- r
}
top_df <- rbindlist(topall)
write.table(top_df, file.path(OUT, "G4_top_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 5. 关键组合并到一张表 -----------------------------------
cat("\n[5] 写汇总\n")
summ <- data.frame(
  signature = names(rankings),
  n_genes   = sapply(rankings, length),
  stringsAsFactors = FALSE
)
write.table(summ, file.path(OUT, "G4_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

cat("\n输出目录：", OUT, "\n")
cat("完成\n")
