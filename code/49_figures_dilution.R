#!/usr/bin/env Rscript
# 49 — Fig11（单细胞样本内配对 + 稀释假说定量排除）与 Fig12（泛癌 IR–微环境地图）
#
# Fig11 回答：bulk 里看到的「IR burden 越高 → 微环境越冷」，
#             会不会只是「免疫细胞 IR 低、浸润多就把 bulk 稀释低了」？
#      A 40 例样本内配对（恶性 vs 免疫）
#      B 三种口径的效应量对比（揭示 Fig6 现行口径被夸大约 10 倍）
#      C 稀释分解：合成 bulk 与浸润比例的关系
#      D 置换检验：纯稀释机制能否造出 TCGA 实测强度
#
# Fig12 回答：这是胶质瘤特有还是泛癌通则？GBM / LGG 处在什么位置？
#      A 31 癌种 × 9 面板的校正后偏 rho 热图
#      B 综合免疫浸润的泛癌排序
#      C 未校正 vs 校正（suppressor 效应在泛癌层面同样成立）
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
OUT <- "output"; OUTP <- "output_pancan"; FIG <- "figures"

S   <- fread(file.path(OUT, "sc_paired_sample.tsv"))
ST  <- fread(file.path(OUT, "sc_paired_stats.tsv"))
DEC <- fread(file.path(OUT, "sc_dilution_decomp.tsv"))
PM  <- fread(file.path(OUT, "sc_dilution_perm.tsv"))
PR  <- fread(file.path(OUT, "sc_dilution_permref.tsv"))
MAP <- fread(file.path(OUTP, "pancan_tme_map.tsv"))

TH <- theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = .2),
        plot.title = element_text(size = 8.3, face = "bold"),
        plot.subtitle = element_text(size = 7.0, colour = "grey35"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7.5),
        legend.position = "none")
CT <- c(Tumor = "#C0392B", Immune = "#185FA5")
GLIO <- c(GBM = "#185FA5", LGG = "#C0392B")

# =============================== Fig11 ===============================
# ---- A：40 例样本内配对斜率图（深度校正后）----
sa <- S[, .(sample, f_immune, Tumor = T_adj, Immune = I_adj)]
sa <- melt(sa, id.vars = c("sample", "f_immune"),
           variable.name = "compartment", value.name = "score")
sa[, compartment := factor(compartment, levels = c("Tumor", "Immune"))]
mm <- sa[, .(m = mean(score)), by = compartment]

pA <- ggplot(sa, aes(compartment, score, group = sample)) +
  geom_line(aes(colour = compartment), linewidth = .18, alpha = .45) +
  geom_point(aes(colour = compartment), size = .7, alpha = .55) +
  geom_line(data = mm, aes(compartment, m, group = 1),
            colour = "grey25", linewidth = .8, inherit.aes = FALSE) +
  geom_point(data = mm, aes(compartment, m), colour = "grey25", size = 2.2,
             inherit.aes = FALSE) +
  scale_colour_manual(values = CT) +
  labs(title = "A  Within-sample pairing (40 tumours)",
       subtitle = sprintf("depth-adjusted; mean difference +%.3f, Wilcoxon p = %.3f",
                          ST[grepl("depth", metric), mean_diff],
                          ST[grepl("depth", metric), wilcoxon_p]),
       x = NULL, y = "IR-axis score (adjusted)") + TH

# ---- B：三种口径的效应量对比（metric 名含空格，必须按内容匹配）----
sb <- ST[grepl("^T", metric)]
sb[grepl("pos", metric),            `:=`(ord = 1, lab = "pos only\n(Fig. 6 current)")]
sb[grepl("depth", metric),          `:=`(ord = 3, lab = "score + depth\nadjusted")]
sb[!grepl("pos|depth", metric),     `:=`(ord = 2, lab = "score\n(internal control)")]
sb[, sig := wilcoxon_p < 0.05]
pB <- ggplot(sb, aes(mean_diff, reorder(lab, ord))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = .3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .12,
                 colour = "grey45", linewidth = .45) +
  geom_point(aes(fill = sig), shape = 21, size = 2.6, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "white")) +
  geom_text(aes(x = ci_hi, label = sprintf("  p = %.2g", wilcoxon_p)),
            size = 2.1, hjust = 0, vjust = .35, colour = "grey25") +
  scale_x_continuous(limits = c(-0.02, 0.085),
                     labels = function(x) sprintf("%.2f", x)) +
  labs(title = "B  Effect size depends on the metric",
       subtitle = "filled = p < 0.05;  whiskers = bootstrap 95% CI",
       x = "mean within-sample difference (tumour - immune)", y = NULL) + TH

# ---- C：稀释分解 ----
S[, bulk := (1 - f_immune) * T_score + f_immune * I_score]
S[, cf_dilute := (1 - f_immune) * mean(T_score) + f_immune * I_score]
S[, cf_intrin := (1 - f_immune) * T_score + f_immune * mean(I_score)]
r_obs <- DEC[grepl("observed|合成", quantity) | quantity == "合成 bulk（观测）", rho_vs_f]
sc <- melt(S[, .(sample, f_immune, `synthetic bulk` = bulk,
                `dilution only` = cf_dilute, `tumour-intrinsic only` = cf_intrin)],
           id.vars = c("sample", "f_immune"),
           variable.name = "scenario", value.name = "score")
pC <- ggplot(sc, aes(f_immune, score)) +
  geom_point(size = .9, alpha = .5, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .5, colour = "grey30") +
  facet_wrap(~ scenario, nrow = 1, scales = "free_y") +
  labs(title = "C  Decomposition of the synthetic bulk signal",
       subtitle = "immune fraction vs synthetic bulk score; lines = linear fit",
       x = "immune cell fraction f", y = NULL) +
  TH + theme(strip.text = element_text(size = 6.4),
             panel.spacing = unit(2.2, "mm"))

# ---- D：置换检验 ----
pd <- PR[, .(reference, tcga_rho, p_dilution_ge)]
pd[, reference := factor(reference, levels = c("LGG IFN-gamma", "LGG CD8_T", "GBM myeloid"))]
pd[, ref_lab := sprintf("%s  (P = %.1f%%)", as.character(reference), p_dilution_ge * 100)]
pD <- ggplot(PM, aes(perm_rho)) +
  geom_histogram(bins = 60, fill = "grey80", colour = "grey55", linewidth = .2) +
  geom_vline(data = pd, aes(xintercept = tcga_rho, colour = ref_lab),
             linewidth = .55) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey40", linewidth = .3) +
  scale_colour_manual(values = c("#C0392B", "#185FA5", "#7F8C8D"),
                      name = "observed TCGA effect (P dilution >=)") +
  labs(title = "D  How strong can pure dilution get?",
       subtitle = "5,000 permutations; T and I independent of f",
       x = expression(paste("Spearman ", rho, " (synthetic bulk vs f)")),
       y = "permutations") +
  TH + theme(legend.position = "top",
             legend.title = element_text(size = 6.2, colour = "grey30"),
             legend.text = element_text(size = 6.2),
             legend.key.width = unit(2.6, "mm"),
             legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 3))

f11 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(3, 5, 3, 3))
ggsave(file.path(FIG, "Fig11_dilution_test.png"), f11,
       width = 200, height = 165, units = "mm", dpi = 600)
ggsave(file.path(FIG, "Fig11_dilution_test.pdf"), f11,
       width = 200, height = 165, units = "mm")

# =============================== Fig12 ===============================
MAPlab <- MAP
MAPlab[, panel_lab := panel]
# 综合免疫放在最左，Sm core（对照）放在最右
pl_order <- c("Immune_aggregate", "T_cell_general", "CD8_T", "Cytotoxic",
              "Myeloid_general", "IFN_gamma_Ayers", "CAF_Stromal",
              "Endothelial", "Sm_core")
MAPlab[, panel_lab := factor(panel_lab, levels = pl_order)]
agg <- MAPlab[panel == "Immune_aggregate"][order(rho_adj)]
MAPlab[, cohort := factor(cohort, levels = agg$cohort)]

pA2 <- ggplot(MAPlab[!is.na(panel_lab)], aes(panel_lab, cohort, fill = rho_adj)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-0.45, 0.45),
                       oob = squish, name = "partial rho") +
  labs(title = "A  IR burden vs microenvironment, 31 tumour types",
       subtitle = "rank partial correlation adjusted for global expression",
       x = NULL, y = NULL) +
  TH + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.2),
             axis.text.y = element_text(size = 5.9),
             legend.position = "right",
             legend.title = element_text(size = 6.5),
             legend.text = element_text(size = 6),
             legend.key.height = unit(8, "mm")) +
  guides(fill = guide_colourbar(barwidth = unit(2.4, "mm")))

agg[, hl := ifelse(cohort %in% c("GBM", "LGG"), cohort, "other")]
agg[, cohort := factor(cohort, levels = rev(cohort))]
pB2 <- ggplot(agg, aes(cohort, rho_adj, fill = hl)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = .3) +
  geom_col(width = .78, colour = NA) +
  scale_fill_manual(values = c(GBM = "#185FA5", LGG = "#C0392B",
                               other = "grey72")) +
  geom_text(aes(y = rho_adj + ifelse(rho_adj < 0, -0.018, 0.018),
                label = ifelse(cohort %in% c("GBM", "LGG"), as.character(cohort), "")),
            size = 2.4, fontface = "bold", colour = "grey20", vjust = .5) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.48, 0.20)) +
  labs(title = "B  Immune-aggregate score, ranked",
       subtitle = sprintf("negative in %d/31 types; q < 0.05 in %d/31",
                          sum(agg$rho_adj < 0), sum(agg$q_adj < 0.05)),
       x = NULL, y = "partial rho (adjusted)") +
  TH + theme(axis.text.y = element_text(size = 5.6),
             panel.grid.major.x = element_line(colour = "grey90", linewidth = .2))

pC2 <- ggplot(MAPlab[panel == "Immune_aggregate"],
              aes(rho_unadj, rho_adj)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey65",
              linewidth = .35) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = .25) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70", linewidth = .25) +
  geom_point(size = 1.5, colour = "grey55") +
  geom_point(data = MAPlab[panel == "Immune_aggregate" &
                             cohort %in% c("GBM", "LGG")],
             aes(colour = cohort), size = 2.8) +
  geom_label(data = MAPlab[panel == "Immune_aggregate" &
                             cohort %in% c("GBM", "LGG")],
             aes(colour = cohort, label = cohort), size = 2.3,
             fill = "white", label.size = NA, label.padding = unit(.6, "mm"),
             vjust = -0.9, fontface = "bold") +
  scale_colour_manual(values = GLIO) +
  annotate("text", x = -0.24, y = 0.14, size = 2.3, colour = "grey25",
           hjust = 0, label = sprintf("%d/31 shift further negative after adjustment\nmedian %.3f -> %.3f",
                                      sum(MAPlab[panel == "Immune_aggregate",
                                                 rho_adj - rho_unadj < 0]),
                                      median(MAPlab[panel == "Immune_aggregate", rho_unadj]),
                                      median(MAPlab[panel == "Immune_aggregate", rho_adj]))) +
  labs(title = "C  Adjustment strengthens, not weakens",
       subtitle = "each point = one tumour type; below the diagonal = stronger",
       x = "unadjusted rho", y = "adjusted rho") + TH

f12 <- (pA2 | (pB2 / pC2)) +
  plot_layout(widths = c(1.15, 1)) &
  theme(plot.margin = margin(3, 5, 3, 3))
ggsave(file.path(FIG, "Fig12_pancan_tme_map.png"), f12,
       width = 205, height = 172, units = "mm", dpi = 600)
ggsave(file.path(FIG, "Fig12_pancan_tme_map.pdf"), f12,
       width = 205, height = 172, units = "mm")

cat("Fig11 / Fig12 written.\n")
