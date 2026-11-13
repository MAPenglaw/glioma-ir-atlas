#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 53 — FigS3：残差化稳健性检验图（评审 P0-1 配套图件）。
#
# 数据源：output_km/residual_robustness.tsv（汇总指标）
#         output_km/figS3_residual_scores.tsv（逐样本 PC1 / burden / TotalExpr）
#
# 面板：
#   A  PC1 方差解释率：R0 原始 vs R1 扣总表达 vs R2 再扣纯度+神经元
#   B  PC1 与总表达量的 Spearman ρ：扣除后归零（证明残差化生效）
#   C  PC1 与 Sm 核心的 ρ：GBM 存活、LGG 仍为零（谱系差异独立于总表达）
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

OUTK <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output_km"
FIG  <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/figures"

d <- fread(file.path(OUTK, "residual_robustness.tsv"))
d[, value := as.numeric(value)]
d <- d[metric %in% c("PC1_var", "rho_PC1_TotalExpr", "rho_PC1_Sm_core",
                     "rho_PC1_burden", "silhouette_k2", "PAC_k2")]
d[, condition := factor(condition,
                        levels = c("R0_raw", "R1_resid_TotalExpr", "R2_resid_full"),
                        labels = c("R0 unadjusted",
                                   "R1 minus global expression",
                                   "R2 also minus purity + neuronal"))]
d[, cancer := factor(cancer, levels = c("GBM", "LGG"))]

COLS <- c("GBM" = "#D85A30", "LGG" = "#378ADD")
TH <- theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "plain", hjust = 0),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        legend.position = "none",
        axis.text.x = element_text(size = 7.5, angle = 18, hjust = 1))

pick <- function(m) d[metric == m]

pA <- ggplot(pick("PC1_var"), aes(condition, value, fill = cancer)) +
  geom_col(position = position_dodge(.75), width = .62,
           colour = "grey25", linewidth = .2) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(.75), vjust = -.5, size = 2.7) +
  scale_fill_manual(values = COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = NULL, y = "PC1 variance explained",
       title = "A  The axis survives removal of global expression",
       subtitle = "R2 retains 84% (GBM) and 85% (LGG) of R0") + TH

pB <- ggplot(pick("rho_PC1_TotalExpr"), aes(condition, value, fill = cancer)) +
  geom_col(position = position_dodge(.75), width = .62,
           colour = "grey25", linewidth = .2) +
  geom_text(aes(label = sprintf("%+.3f", value)),
            position = position_dodge(.75), vjust = -.5, size = 2.7) +
  geom_hline(yintercept = 0, linewidth = .35, colour = "grey40") +
  scale_fill_manual(values = COLS) +
  scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
  labs(x = NULL, y = expression(rho ~ "(PC1, global expression)"),
       title = "B  Residualisation works: the association is abolished",
       subtitle = "R0 +0.558 / +0.348  ->  R2 +0.030 / +0.012") + TH

pC <- ggplot(pick("rho_PC1_Sm_core"), aes(condition, value, fill = cancer)) +
  geom_col(position = position_dodge(.75), width = .62,
           colour = "grey25", linewidth = .2) +
  geom_text(aes(label = sprintf("%+.3f", value)),
            position = position_dodge(.75), vjust = -.5, size = 2.7) +
  geom_hline(yintercept = 0, linewidth = .35, colour = "grey40") +
  scale_fill_manual(values = COLS) +
  scale_y_continuous(expand = expansion(mult = c(.15, .22))) +
  labs(x = NULL, y = expression(rho ~ "(PC1, Sm core)"),
       title = "C  The spliceosome coupling survives in GBM only",
       subtitle = "GBM -0.258 (p = 1.3e-3); LGG +0.025 (p = 0.58)") + TH

p <- (pA | pB | pC) +
  plot_annotation(
    title = "Robustness of the continuous IR axis to removal of global transcript output",
    subtitle = "Every high-variability IR event residualised on the listed covariates (rank-normalised); all diagnostics recomputed",
    theme = theme(plot.title = element_text(size = 11.5),
                  plot.subtitle = element_text(size = 9, colour = "grey35")))

ggsave(file.path(FIG, "FigS3_residual_robustness.pdf"), p, width = 12.4, height = 4.5)
ggsave(file.path(FIG, "FigS3_residual_robustness.png"), p, width = 12.4, height = 4.5,
       dpi = 300)
cat("[sav] FigS3_residual_robustness\n[done]\n")
