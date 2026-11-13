#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 36 — 绘制泛癌验证图（Fig8）与分子分型阴性图（FigS1）。
#
# Fig8  泛癌验证：31 个癌种上「IR 是连续轴」的四项证据
# FigS1 IR 轴独立于 IDH / 1p19q / WHO 分级（阴性对照）
#
# 数据源：output_pancan/pancan_diagnostics.tsv、k2_bisection.tsv、
#         mol_burden_plot.tsv、mol_burden_stats.tsv
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output_pancan"
FIG <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/figures"

diag <- fread(file.path(OUT, "pancan_diagnostics.tsv"))
bis  <- fread(file.path(OUT, "k2_bisection.tsv"))
dt   <- merge(diag[, .(projectID, n_samples, PC1_evr, rho_PC1_burden,
                       silhouette_k2, PAC_k2, BIC_prefers_k1)],
              bis[, .(projectID, AUC_k2label_vs_burden,
                      ARI_k2label_vs_median_split)],
              by = "projectID", all = TRUE)

dt[, short := sub("TCGA-", "", projectID)]
dt[, glio := fifelse(short %in% c("GBM", "LGG"), short, "other")]
dt$glio <- factor(dt$glio, levels = c("GBM", "LGG", "other"))

# 按样本量排序，便于阅读
dt <- dt[order(-n_samples)]
dt$short <- factor(dt$short, levels = dt$short)

COLS <- c("GBM" = "#D85A30", "LGG" = "#378ADD", "other" = "#888780")

TH <- theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 12, face = "plain", hjust = 0),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "none")

dotp <- function(d, xvar, xlab, title, sub, vline = NA, vlab = "") {
  p <- ggplot(d, aes(x = .data[[xvar]], y = short, colour = glio)) +
    geom_vline(xintercept = vline, linetype = 2, colour = "grey60",
               linewidth = .3) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = COLS) +
    scale_x_continuous(labels = number_format(accuracy = .01)) +
    labs(x = xlab, y = NULL, title = title, subtitle = sub) + TH
  if (!is.na(vline)) {
    p <- p + annotate("text", x = vline, y = 1, label = vlab,
                      size = 2.6, colour = "grey40", vjust = -0.4, hjust = 0)
  }
  p
}

pA <- dotp(dt, "PAC_k2",
           "PAC (k=2, corrected)", "A  Ambiguous-clustering proportion",
           "Median 0.535; read with silhouette — low PAC alone does not imply subtypes")

pB <- dotp(dt, "rho_PC1_burden",
           "Spearman rho(PC1, IR burden)",
           "B  PC1 is the burden axis itself",
           "Median 0.95 — the leading axis adds nothing beyond burden")

pC <- dotp(dt, "AUC_k2label_vs_burden",
           "AUC (k=2 label discriminates burden)",
           "C  The two 'clusters' are just the axis cut in half",
           "27/31 tumour types >= 0.9 (k=2 bisects a continuum)",
           vline = 0.9, vlab = " 0.9")

pD <- ggplot(dt, aes(x = silhouette_k2, y = short, colour = glio)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60", linewidth = .3) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = COLS) +
  scale_x_continuous(limits = c(0, 0.62), expand = c(0.02, 0.02)) +
  labs(x = "Silhouette (k=2)", y = NULL,
       title = "D  Cluster separation is weak everywhere",
       subtitle = "Dashed line 0.5 = strong structure. Only 1/31 (KICH, 0.508) marginally exceeds it; median 0.28") + TH

f8 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "A continuous intron-retention axis is a pan-cancer property",
    subtitle = sprintf(
      "%d TCGA tumour types, %s primary tumours, %s intron-retention events with >=90%% pan-cancer detection.\nGBM (n=%d) and LGG (n=%d) highlighted; grey = other tumour types.",
      nrow(dt), format(sum(dt$n_samples), big.mark = ","), "34,752",
      dt[short == "GBM", n_samples], dt[short == "LGG", n_samples]),
    theme = theme(plot.title = element_text(size = 13, face = "plain"),
                  plot.subtitle = element_text(size = 9, colour = "grey35")))

ggsave(file.path(FIG, "Fig8_pancan_continuous_axis.png"), f8,
       width = 12.5, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(FIG, "Fig8_pancan_continuous_axis.pdf"), f8,
       width = 12.5, height = 9.5, bg = "white")
cat("[sav] Fig8_pancan_continuous_axis\n")

# ------------------------------------------------------------------ FigS1
mb <- fread(file.path(OUT, "mol_burden_plot.tsv"))
st <- fread(file.path(OUT, "mol_burden_stats.tsv"))
mb <- mb[!is.na(IDH_status)]

p1 <- ggplot(mb[cohort == "LGG"], aes(IDH_status, RI_burden, fill = IDH_status)) +
  geom_boxplot(width = .45, outlier.shape = NA, alpha = .8) +
  geom_jitter(width = .12, size = .8, alpha = .45) +
  scale_fill_manual(values = c("IDHmut" = "#1D9E75", "IDHwt" = "#D85A30")) +
  labs(title = "A  LGG: IR burden by IDH status",
       subtitle = "Cliff's delta 0.067 (95% CI -0.09 to 0.23), p = 0.40",
       x = NULL, y = "IR burden (mean PSI)") + TH

p2 <- ggplot(mb[cohort == "LGG" & !is.na(mol3)], aes(mol3, RI_burden, fill = mol3)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .8) +
  geom_jitter(width = .12, size = .8, alpha = .45) +
  scale_fill_manual(values = c("codel_like" = "#378ADD",
                               "noncodel_like" = "#1D9E75",
                               "IDHwt" = "#D85A30")) +
  labs(title = "B  LGG: three molecular classes",
       subtitle = "Kruskal-Wallis p = 0.70 (burden) / 0.97 (PC1)",
       x = NULL, y = "IR burden (mean PSI)") + TH +
  theme(axis.text.x = element_text(size = 8, angle = 15, hjust = 1))

p3 <- ggplot(mb[cohort == "GBM"], aes(IDH_status, RI_burden, fill = IDH_status)) +
  geom_boxplot(width = .45, outlier.shape = NA, alpha = .8) +
  geom_jitter(width = .12, size = 1.1, alpha = .7) +
  scale_fill_manual(values = c("IDHmut" = "#1D9E75", "IDHwt" = "#D85A30")) +
  labs(title = "C  GBM: IR burden by IDH status",
       subtitle = "Only 7 IDH-mutant GBM — underpowered, delta 0.05 [-0.46, 0.57]",
       x = NULL, y = "IR burden (mean PSI)") + TH

fr <- st[!is.na(cliffs_delta) & comparison != "Kruskal 3-class (codel/noncodel/IDHwt)"]
fr <- fr[, metric := ifelse(metric == "RI_burden", "IR burden", metric)]
fr <- fr[, term := paste0(cohort, ": ", comparison, " (", metric, ")")]
fr <- fr[order(term)]
fr$term <- factor(fr$term, levels = rev(fr$term))

pfmt <- function(p) ifelse(is.na(p), "NA",
                    ifelse(p < 1e-4, sprintf("p = %.1e", p),
                           sprintf("p = %.2f", p)))

p4 <- ggplot(fr, aes(term, cliffs_delta)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .12, linewidth = .35) +
  geom_point(size = 2.4, colour = "#534AB7") +
  coord_flip(ylim = c(-0.7, 1.18)) +
  geom_text(aes(y = ci_hi + .05,
                label = sprintf("%.3f (%s)", cliffs_delta, pfmt(p))),
            size = 2.7, hjust = 0) +
  labs(x = NULL, y = "Cliff's delta (95% CI)",
       title = "D  Effect sizes all straddle zero",
       subtitle = "No comparison excludes even a small effect in either direction") +
  TH + theme(axis.text.y = element_text(size = 7.5))

fs1 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "The IR axis is independent of IDH, 1p/19q and WHO grade",
    subtitle = "575 TCGA gliomas with both IR burden and molecular annotation (GBM n=110, LGG n=465). IDH status from IDH1/IDH2 mutation; 1p/19q codeletion proxied by CIC (19q13.2) / FUBP1 (1p31.1) mutation.",
    theme = theme(plot.title = element_text(size = 13, face = "plain"),
                  plot.subtitle = element_text(size = 9, colour = "grey35")))

ggsave(file.path(FIG, "FigS1_molecular_independence.png"), fs1,
       width = 12.5, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(FIG, "FigS1_molecular_independence.pdf"), fs1,
       width = 12.5, height = 8.5, bg = "white")
cat("[sav] FigS1_molecular_independence\n")
cat("[done]\n")
