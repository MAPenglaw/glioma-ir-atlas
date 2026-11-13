#!/usr/bin/env Rscript
# Step 7 — 主图（发表级，英文标签）
suppressMessages({
  library(ggplot2); library(data.table); library(patchwork)
  library(survival); library(survminer); library(reshape2)
})

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
FIG <- paste0(OUT, "/figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
TH <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(size = 11, face = "bold"))

tab <- read.delim(paste0(OUT, "/output/sample_RI_group_annotated.tsv"),
                  row.names = 1, check.names = FALSE)
ann <- read.delim(paste0(OUT, "/output/sample_RI_group_clinical.tsv"),
                  row.names = 1, check.names = FALSE)
tab$Subtype <- c("0"="Classical","1"="Neural","2"="Proneural","3"="Mesenchymal")[
  as.character(ann[rownames(tab), "GeneExp_Subtype"])]
tab$RI_group <- factor(tab$RI_group, levels = c("RI-low", "RI-high"))
pal <- c("RI-low" = "#4C72B0", "RI-high" = "#C44E52")

## =================== Fig 1: RI 图谱总览 ===================
p1a <- ggplot(tab, aes(x = RI_burden)) +
  geom_histogram(bins = 24, fill = "#4C72B0", colour = "white", linewidth = 0.2) +
  geom_vline(aes(xintercept = mean(RI_burden)), linetype = 2, colour = "grey30") +
  labs(x = "Global RI burden (mean PSI of detected IR events)", y = "Samples",
       title = "A  Continuous distribution of RI burden",
       subtitle = sprintf("n = %d TCGA-GBM; range %.3f-%.3f", nrow(tab),
                          min(tab$RI_burden), max(tab$RI_burden))) + TH

p1b <- ggplot(tab, aes(x = RI_PC1, y = RI_PC2, colour = RI_group)) +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_manual(values = pal) +
  labs(x = "PC1 (25.3% variance)", y = "PC2 (3.7%)",
       title = "B  PCA of 1,672 high-variance IR events",
       subtitle = "single dominant axis, no separated islands") + TH +
  theme(legend.position = c(0.82, 0.82), legend.background = element_rect(fill = "white", colour = NA))

kd <- c(2, 3, 4, 5, 6)
grid <- seq(0, 1, length.out = 101)
cdf_dt <- rbindlist(lapply(kd, function(k) {
  C <- as.matrix(read.delim(paste0(OUT, "/output/consensus_k", k, ".tsv"), row.names = 1, check.names = FALSE))
  v <- C[upper.tri(C)]
  data.table(k = factor(k), x = grid,
             y = sapply(grid, function(g) mean(v <= g)))
}))
p1c <- ggplot(cdf_dt, aes(x, y, colour = k)) + geom_line(linewidth = 0.7) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Consensus index", y = "CDF", colour = "k",
       title = "C  Consensus CDF (no clear elbow)",
       subtitle = "PAC 0.52 at k=2, no elbow; max silhouette 0.14") + TH +
  theme(legend.position = c(0.85, 0.3), legend.key.size = unit(9, "pt"))

cons <- as.matrix(read.delim(paste0(OUT, "/output/consensus_k2.tsv"), row.names = 1, check.names = FALSE))
ord <- order(tab$RI_burden)
hm <- as.data.frame(as.table(cons[ord, ord])); colnames(hm) <- c("Var1","Var2","value")
p1d <- ggplot(hm, aes(Var1, Var2, fill = value)) + geom_raster() +
  scale_fill_gradient(low = "white", high = "#C44E52", limits = c(0, 1), name = "consensus") +
  coord_fixed() + labs(x = NULL, y = NULL,
       title = "D  Consensus matrix, k = 2",
       subtitle = "samples ordered by RI burden; no block structure") + TH +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

f1 <- (p1a | p1b) / (p1c | p1d)
ggsave(paste0(FIG, "/Fig1_RI_landscape.png"), f1, width = 9.5, height = 8, dpi = 300)
pdf(paste0(FIG, "/Fig1_RI_landscape.pdf"), width = 9.5, height = 8, onefile = TRUE)
print(f1); dev.off()
cat("Fig1 done\n")

## =================== Fig 2: RI 轴的生物学属性 ===================
co <- read.delim(paste0(OUT, "/output/step3a_confound_corr.tsv"))
co$variable <- factor(co$variable, levels = co$variable[order(co$rho_burden)])
p2a <- ggplot(co, aes(rho_burden, variable, fill = abs(rho_burden) > 0.4)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "#9FB6CD", "TRUE" = "#C44E52"), guide = "none") +
  labs(x = "Spearman rho with RI burden", y = NULL,
       title = "A  RI burden vs potential confounders",
       subtitle = "proliferation / purity / stromal all |rho| < 0.05") + TH

g <- read.delim(paste0(OUT, "/output/gsva_scores.tsv"), row.names = 1, check.names = FALSE)
dd <- data.frame(burden = tab$RI_burden,
                 U1 = g[rownames(tab), "U1_snRNP"],
                 SPL = g[rownames(tab), "Spliceosome_core"])
p2b <- ggplot(dd, aes(burden, U1)) + geom_point(size = 1.6, colour = "#4C72B0", alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#C44E52", linewidth = .7) +
  labs(x = "RI burden", y = "U1 snRNP score (ssGSEA)",
       title = "B  Higher RI tracks loss of U1 snRNP",
       subtitle = sprintf("rho = %.2f, BH p = 5.5e-03",
         cor(dd$burden, dd$U1, method = "spearman"))) + TH
p2c <- ggplot(dd, aes(burden, SPL)) + geom_point(size = 1.6, colour = "#55A868", alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#C44E52", linewidth = .7) +
  labs(x = "RI burden", y = "Core spliceosome score (ssGSEA)",
       title = "C  ... with a reciprocal rise in core spliceosome",
       subtitle = sprintf("rho = %.2f, BH p = 6.2e-03",
         cor(dd$burden, dd$SPL, method = "spearman"))) + TH

gc <- read.delim(paste0(OUT, "/output/gsva_group_compare.tsv"))
gc <- gc[order(gc$rho_burden), ]
gc$gene_set <- factor(gc$gene_set, levels = gc$gene_set)
gc$sig <- ifelse(gc$p_adj < 0.05, "BH < 0.05", "n.s.")
p2d <- ggplot(gc, aes(rho_burden, gene_set, colour = sig)) +
  geom_segment(aes(x = 0, xend = rho_burden, yend = gene_set), linewidth = .5) +
  geom_point(size = 2) + geom_vline(xintercept = 0, linewidth = .3) +
  scale_colour_manual(values = c("BH < 0.05" = "#C44E52", "n.s." = "#9E9E9E")) +
  labs(x = "Spearman rho with RI burden", y = NULL, colour = NULL,
       title = "D  Immune programmes are largely unchanged",
       subtitle = "only U1 snRNP + spliceosome core survive BH") + TH +
  theme(legend.position = c(0.83, 0.16),
        legend.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 10.5, face = "bold"),
        plot.subtitle = element_text(size = 9))

f2 <- (p2a | p2b) / (p2c | p2d)
ggsave(paste0(FIG, "/Fig2_RI_axis_biology.png"), f2, width = 10, height = 9, dpi = 300)
pdf(paste0(FIG, "/Fig2_RI_axis_biology.pdf"), width = 10, height = 9, onefile = TRUE)
print(f2); dev.off()
cat("Fig2 done\n")

## =================== Fig 3: 预后价值（阴性结果）===================
sv <- readRDS(paste0(OUT, "/output/survival_input.rds"))
sv$RI_group <- factor(sv$RI_group, levels = c("RI-low", "RI-high"))
fit_os <- survfit(Surv(OS.time, OS) ~ RI_group, data = sv)
fit_pfi <- survfit(Surv(PFI.time, PFI) ~ RI_group, data = sv)
km <- function(f, lab, sub) {
  ggsurvplot(f, data = sv, pval = TRUE, pval.size = 3.2, conf.int = FALSE,
             palette = unname(pal), risk.table = TRUE, tables.height = 0.26,
             legend.title = "", legend.labs = c("RI-low", "RI-high"),
             xlab = "Time (days)", ylab = lab, title = sub,
             ggtheme = TH, size = 0.8, censor.size = 2)
}
k1 <- km(fit_os, "Overall survival", "A  OS by RI group (n = 150)")$plot
k2 <- km(fit_pfi, "Progression-free interval", "B  PFI by RI group")$plot

bat <- read.delim(paste0(OUT, "/output/survival_metric_battery.tsv"))
bat <- bat[bat$endpoint == "OS", ]
bat$metric <- factor(bat$metric, levels = bat$metric[order(bat$HR_per_SD)])
bat$sig <- ifelse(bat$p_adj_within_endpoint < 0.05, "BH < 0.05", "n.s.")
p3c <- ggplot(bat, aes(HR_per_SD, metric, colour = sig)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .25, linewidth = .5) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("BH < 0.05" = "#C44E52", "n.s." = "#4C72B0")) +
  coord_cartesian(xlim = c(0.6, 1.45)) +
  labs(x = "Hazard ratio per +1 SD (OS)", y = NULL, colour = NULL,
       title = "C  16 RI-derived metrics, none prognostic",
       subtitle = "positive control: age HR = 1.027/yr, p = 6e-04") + TH +
  theme(legend.position = c(0.82, 0.12))

ca <- readRDS(paste0(OUT, "/output/cgga_analysis.rds"))
cbu <- read.delim(paste0(OUT, "/output/cgga_cox_univariate.tsv"))
cbu <- cbu[!is.na(cbu$subset) & nzchar(trimws(cbu$subset)), ]
cbu <- cbu[!cbu$subset %in% c("1p19q codel", "1p19q non-codel"), ]
cbu$subset <- factor(cbu$subset, levels = cbu$subset[order(cbu$HR_per_SD)])
p3d <- ggplot(cbu, aes(HR_per_SD, subset)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .25, linewidth = .5, colour = "#4C72B0") +
  geom_point(size = 2, colour = "#4C72B0") +
  labs(x = "Hazard ratio per +1 SD (CGGA, n = 970)", y = NULL,
       title = "D  CGGA: signal confined to IDH-mutant / low grade",
       subtitle = "attenuated after adjusting grade+IDH (HR 0.97, p=0.54)") + TH

f3 <- (k1 | k2) / (p3c | p3d)
ggsave(paste0(FIG, "/Fig3_prognostic_null.png"), f3, width = 11, height = 9, dpi = 300)
pdf(paste0(FIG, "/Fig3_prognostic_null.pdf"), width = 11, height = 9, onefile = TRUE)
print(f3); dev.off()
cat("Fig3 done\n")

## =================== Fig 4: 机制锚定 GSE227187 ===================
ev <- data.frame(type = c("SE", "A5SS", "A3SS", "MXE", "RI"),
                 up = c(2485, 802, 940, 1024, 794),
                 down = c(2764, 670, 776, 1176, 606))
evl <- reshape2::melt(ev, id.vars = "type", variable.name = "dir", value.name = "n")
evl$dir <- factor(evl$dir, levels = c("down", "up"),
                  labels = c("inclusion down", "inclusion up"))
p4a <- ggplot(evl, aes(type, n, fill = dir)) + geom_col(position = "dodge", width = .72) +
  scale_fill_manual(values = c("inclusion down" = "#9FB6CD", "inclusion up" = "#C44E52")) +
  labs(x = NULL, y = "Significant events (FDR<0.05, |dPSI|>0.1)", fill = NULL,
       title = "A  SNRPC knockdown reshuffles the splicing programme",
       subtitle = "GSE227187, sgSNRPC vs sgControl (2v2)") + TH

ox <- data.frame(
  set = factor(c("OXPHOS all", "Complex I", "Complex III", "Complex IV", "Complex V",
                 "Mito ribosome"),
               levels = rev(c("OXPHOS all", "Complex I", "Complex III", "Complex IV",
                              "Complex V", "Mito ribosome"))),
  n = c(167, 140, 16, 51, 24, 106),
  med = c(-0.157, -0.148, -0.162, -0.148, -0.259, -0.111),
  pct = c(77, 77, 94, 75, 79, 72),
  p = c(5.8e-17, 5.3e-14, 2.1e-04, 8.6e-05, 1.5e-04, 7.7e-11))
p4b <- ggplot(ox, aes(med, set)) +
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_col(aes(fill = pct), width = .68) +
  geom_text(aes(label = sprintf("%d%% down, n=%d", pct, n)), hjust = -0.08, size = 2.9) +
  scale_fill_gradient(low = "#F2C1C1", high = "#A32D2D", name = "% down") +
  coord_cartesian(xlim = c(-0.32, 0.06)) +
  labs(x = "Median log2 fold change (SNRPC KD)", y = NULL,
       title = "B  OXPHOS genes are coordinately down-regulated",
       subtitle = "167 OXPHOS genes, 77% down, Wilcoxon p = 5.8e-17") + TH +
  theme(legend.position = "none")

nmd <- data.frame(g = c("Expression down\n(NMD-consistent)", "Expression up",
                        "RI retention up, no change"),
                  n = c(21, 9, 0))
p4c <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  annotate("rect", xmin = .5, xmax = 5.5, ymin = .55, ymax = .95,
           fill = "#EAF3DE", colour = "#639922", linewidth = .4) +
  annotate("rect", xmin = .5, xmax = 5.5, ymin = .1, ymax = .5,
           fill = "#F1EFE8", colour = "#B4B2A9", linewidth = .4) +
  annotate("text", x = 3, y = .82, size = 3.4, fontface = "bold", colour = "#27500A",
           label = "21 / 30 RI-retention genes (70%) show reduced mRNA") +
  annotate("text", x = 3, y = .66, size = 3,
           label = "supports intron retention -> NMD -> down-regulation") +
  annotate("text", x = 3, y = .30, size = 3.1, colour = "#444441",
           label = "Representative: UQCRC1 RI +0.22 / log2FC -0.37") +
  annotate("text", x = 3, y = .16, size = 3.1, colour = "#444441",
           label = "COX7C RI +0.33 / log2FC -0.15;  MRS2 RI +0.65") +
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 1)) + theme_void() +
  labs(title = "C  Intron retention is coupled to transcript loss") +
  theme(plot.title = element_text(size = 11, face = "bold"))

cta <- read.delim(paste0(OUT, "/output/sc_celltype_RIaxis.tsv"), row.names = 1)
cta$cell_type <- factor(rownames(cta), levels = rownames(cta)[order(cta$RIaxis_pos)])
p4d <- ggplot(cta, aes(RIaxis_pos, cell_type, fill = cell_type)) +
  geom_col(width = .7) + scale_fill_brewer(palette = "Set2", guide = "none") +
  geom_text(aes(label = sprintf("n = %s", format(n_cells, big.mark = ","))),
            hjust = -0.05, size = 2.8) +
  coord_cartesian(xlim = c(0, max(cta$RIaxis_pos) * 1.45)) +
  labs(x = "Mean expression of RI-axis positive genes (log1p CP10K)", y = NULL,
       title = "D  Single-cell localization (GSE182109, 264,523 cells)",
       subtitle = "RI-axis positive genes peak in malignant cells (Mann-Whitney p < 1e-300)") + TH

f4 <- (p4a | p4b) / (p4c | p4d)
ggsave(paste0(FIG, "/Fig4_mechanism_anchor.png"), f4, width = 11, height = 8.5, dpi = 300)
pdf(paste0(FIG, "/Fig4_mechanism_anchor.pdf"), width = 11, height = 8.5, onefile = TRUE)
print(f4); dev.off()
cat("Fig4 done\n")

cat("\n全部图完成 -> ", FIG, "\n")
