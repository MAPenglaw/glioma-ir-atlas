#!/usr/bin/env Rscript
# Step 3b — GSVA(ssGSEA) 免疫微环境与干性分析
# 基因集: 自建 marker panel（全部显式列出，便于复现与写进方法学）
# 方法: GSVA::gsva(method="ssgsea", kcdf="Gaussian")  Barbie 2009 / Hänzelmann 2013

suppressMessages({
  library(GSVA); library(data.table)
})

SRC <- paste0("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/",
              "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
bc <- function(x) substr(gsub('"', '', x), 1, 15)

PANELS <- list(
  CD8_T_cell          = c("CD8A","CD8B","GZMA","GZMK","PRF1","CD3D","CD3E","CD2","IL2RB","NKG7"),
  CD4_T_cell          = c("CD4","CD3D","CD3E","CD40LG","IL7R","CCR7","LTB","TRAT1","ICOS","CD28"),
  Treg                = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TNFRSF4","CCR8","BATF","LAYN","IL10"),
  NK_cell             = c("KLRD1","KLRF1","NKG7","GNLY","GZMB","NCAM1","NCR1","PRF1","KLRB1","FCGR3A"),
  B_cell              = c("CD19","MS4A1","CD79A","CD79B","PAX5","TCL1A","BLNK","CD27","IGKC","TNFRSF17"),
  Monocyte            = c("S100A8","S100A9","S100A12","CD14","LYZ","VCAN","FCN1","CST3","CLEC7A","FCGR3A"),
  Macrophage_M1       = c("NOS2","IL1B","TNF","CXCL9","CXCL10","CXCL11","SOCS3","CD80","CD86","IDO1"),
  Macrophage_M2       = c("CD163","MRC1","MSR1","CD209","MAF","IGF1","CCL13","CCL18","TGFB1","IL10"),
  Neutrophil          = c("FCGR3B","CSF3R","CXCR2","ELANE","MPO","DEFA3","CEACAM8","FPR1","S100A12","MMP9"),
  Dendritic_cell      = c("CD1C","CD1E","CLEC9A","CLEC10A","BATF3","ITGAX","FCER1A","LILRA4","IRF8","LAMP3"),
  Mast_cell           = c("TPSAB1","TPSB2","CPA3","MS4A2","KIT","CMA1","HDC","GATA2","SIGLEC6","CTSG"),
  Immune_checkpoint   = c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","VSIR","IDO1","CD276","BTLA"),
  IFN_gamma_signature = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG","ISG20","GBP1","CD274","CXCR6"),
  Cytolytic_activity  = c("GZMA","PRF1"),
  Tcell_inflamed_TIS  = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
                          "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
                          "PSMB10","STAT1","TIGIT"),
  Stromal_fibrosis    = c("FAP","POSTN","COL1A1","COL1A2","COL3A1","FN1","DCN","LUM","THBS2",
                          "ACTA2","PDGFRB","TIMP1"),
  Proliferation       = c("MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","BUB1","BUB1B","AURKA",
                          "AURKB","PLK1","CENPA","CENPF","RRM2","TYMS","MCM2","MCM6","E2F1","FOXM1","TTK"),
  Glioma_StemCell     = c("SOX2","NANOG","POU5F1","PROM1","NES","ALDH1A3","OLIG2","SALL2",
                          "BMI1","KLF4","MYC","FUT4","ITGA6","THY1","PTPRZ1"),
  U1_snRNP            = c("SNRPC","SNRPA","SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF",
                          "SNRPG","SNRNP70","SNRNP200","PRPF8"),
  Spliceosome_core    = c("SF3B1","SF3B2","SF3A1","U2AF1","U2AF2","SRSF1","SRSF2","SRSF3",
                          "HNRNPA1","HNRNPC","PTBP1","RBM39")
)

## ---------- 载入表达谱 ----------
expr <- fread(paste0(SRC, "/GBM_HiSeqV2.gz"), sep = "\t", header = TRUE,
              data.table = FALSE, check.names = FALSE)
rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
colnames(expr) <- bc(colnames(expr))
expr <- expr[, !duplicated(colnames(expr))]
expr <- as.matrix(expr)

tab <- read.delim(paste0(OUT, "/output/sample_RI_group_annotated.tsv"),
                  row.names = 1, check.names = FALSE)
samples <- intersect(rownames(tab), colnames(expr))
expr <- expr[, samples]
tab <- tab[samples, ]
cat("表达矩阵:", nrow(expr), "基因 ×", ncol(expr), "样本\n")

## ---------- GSVA / ssGSEA ----------
geneSets <- lapply(PANELS, function(g) intersect(g, rownames(expr)))
nG <- sapply(geneSets, length)
cat("基因集命中情况:\n"); print(nG)
geneSets <- geneSets[nG >= 3]

set.seed(20260917)
## GSVA >= 2.x 采用 Param 对象接口：ssgseaParam 即标准 ssGSEA（Barbie 2009）
res <- gsva(ssgseaParam(expr, geneSets, alpha = 0.25, normalize = TRUE,
                        verbose = FALSE))
cat("GSVA 输出:", nrow(res), "基因集 ×", ncol(res), "样本  (class:", class(res)[1], ")\n")
res <- as.matrix(res)
storage.mode(res) <- "double"
colnames(res) <- samples
write.table(t(res), paste0(OUT, "/output/gsva_scores.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)

## ---------- 组间比较 ----------
sc <- read.delim(paste0(OUT, "/output/gsva_scores.tsv"), row.names = 1,
                 check.names = FALSE)
sc[] <- lapply(sc, as.numeric)
sc$RI_group <- as.character(tab[rownames(sc), "RI_group"])
sc$RI_burden <- as.numeric(tab[rownames(sc), "RI_burden"])
GENE_SETS <- setdiff(colnames(sc), c("RI_group", "RI_burden"))

cat("\n================ RI-high vs RI-low (Wilcoxon) ================\n")
cat(sprintf("%-22s %9s %9s %11s %11s\n", "gene_set", "high", "low", "diff", "p"))
out <- data.frame()
for (g in GENE_SETS) {
  a <- sc[sc$RI_group == "RI-high", g]
  b <- sc[sc$RI_group == "RI-low",  g]
  w <- suppressWarnings(wilcox.test(a, b))
  r <- suppressWarnings(cor.test(sc$RI_burden, sc[[g]], method = "spearman"))
  cat(sprintf("%-22s %9.3f %9.3f %+11.3f %11.2e\n",
              g, mean(a), mean(b), mean(a) - mean(b), w$p.value))
  out <- rbind(out, data.frame(gene_set = g, mean_RIhigh = mean(a),
                               mean_RIlow = mean(b), diff = mean(a) - mean(b),
                               p_wilcox = w$p.value,
                               rho_burden = as.numeric(r$estimate),
                               p_rho = r$p.value))
}
out$p_adj <- p.adjust(out$p_wilcox, method = "BH")
out <- out[order(out$p_wilcox), ]
write.table(out, paste0(OUT, "/output/gsva_group_compare.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nBH 校正后 p<0.05 的基因集数:", sum(out$p_adj < 0.05), "/", nrow(out), "\n")
print(out[out$p_adj < 0.05, c("gene_set","diff","p_wilcox","p_adj","rho_burden")],
      row.names = FALSE)

cat("\n完成 -> gsva_scores.tsv / gsva_group_compare.tsv\n")
