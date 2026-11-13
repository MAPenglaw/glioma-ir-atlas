#!/usr/bin/env Rscript
# 46 — Fig10：肿瘤微环境耦合的混杂排查（偏相关 + 总表达量诊断）
# 回答：LGG 的「IR burden 越高 → 免疫越冷」是不是总表达量 / 神经元 /
#       剪接体 / 纯度 / IDH 的影子？GBM 为什么在未校正时看不见免疫信号？
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
OUT <- "output_pancan"; FIG <- "figures"
pc <- fread(file.path(OUT, "tme_partial_coupling.tsv"))
bh <- fread(file.path(OUT, "tme_partial_bh.tsv"))
cf <- fread(file.path(OUT, "tme_confounders.tsv"))
sm <- unique(fread(file.path(OUT, "tme_partial_samples.tsv")),
             by = c("sample", "cohort"))

TH <- theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = .2),
        plot.title = element_text(size = 8.3, face = "bold"),
        plot.subtitle = element_text(size = 7.3, colour = "grey35"),
        axis.text = element_text(size = 7),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 7.5),
        legend.position = "none")

# 只看微环境面板（技术/对照面板在 C 面板里单独展示）
TME <- c("T_cell_general", "CD8_T", "Treg", "NK", "B_cell", "Myeloid_general",
         "Microglia", "M2_like", "Neutrophil", "DC", "CAF_Stromal",
         "Endothelial", "IFN_gamma_Ayers", "Cytotoxic", "MHC_I", "MHC_II",
         "Checkpoint")
M0 <- "M0_raw"; M7 <- "M7_hi_lo_layers"
COH_COL <- c(GBM = "#185FA5", LGG = "#C0392B")

# ---------- A / B：哑铃图（未校正 -> 最保守校正）----------
panel_dumbbell <- function(coh, ttl, sub) {
  s <- pc[cohort == coh & panel %in% TME & model %in% c(M0, M7)]
  s <- merge(s, bh[cohort == coh, .(panel, model, q)], by = c("panel", "model"))
  s <- dcast(s, panel ~ model, value.var = c("rho", "p", "q"))
  setorderv(s, "rho_M7_hi_lo_layers")
  s[, panel := factor(panel, levels = panel)]
  # 填充 = BH 校正后 q < 0.05（与图注一致）
  s[, sig7 := q_M7_hi_lo_layers < 0.05]
  ggplot(s, aes(y = panel)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = .3) +
    geom_segment(aes(x = rho_M0_raw, xend = rho_M7_hi_lo_layers,
                     yend = panel), colour = "grey75", linewidth = .35) +
    geom_point(aes(x = rho_M0_raw), shape = 1, size = 1.9,
               colour = "grey55") +
    geom_point(aes(x = rho_M7_hi_lo_layers, fill = sig7),
               shape = 21, size = 2.3, colour = "black", stroke = .25) +
    scale_fill_manual(values = c(`TRUE` = COH_COL[[coh]], `FALSE` = "white")) +
    scale_x_continuous(expand = expansion(mult = c(.10, .16)),
                       limits = c(-0.52, 0.22)) +
    labs(x = "Spearman rho with IR burden", y = NULL,
         title = ttl, subtitle = sub) + TH
}

nA <- bh[cohort == "LGG" & model == M7 & panel %in% TME & q < 0.05, .N]
nB <- bh[cohort == "GBM" & model == M7 & panel %in% TME & q < 0.05, .N]
pA <- panel_dumbbell("LGG", "A  LGG (n=498): survives adjustment",
                     sprintf("open = unadjusted; filled = M7; red = q<0.05 (%d/17, all negative)", nA))
pB <- panel_dumbbell("GBM", "B  GBM (n=150): emerges after adjustment",
                     sprintf("0/17 unadjusted (MDE 0.23); %d/17 at q<0.05", nB))

# ---------- C：总表达量混杂有多强 ----------
cf[, cohort := factor(cohort, levels = c("GBM", "LGG"))]
rho_lab <- pc[model == M0 & panel %in% c("TotalExpr", "TotalExpr_hi", "TotalExpr_lo"),
              .(cohort, panel,
                lab = sprintf("%s ~ burden: rho=%+.2f",
                              fifelse(panel == "TotalExpr", "all",
                                      fifelse(panel == "TotalExpr_hi", "high-expr",
                                              "low-expr")), rho))]
lab_all <- rho_lab[panel == "TotalExpr"][, .(cohort, lab = sub("all ~ burden: ", "", lab))]
cf[lab_all, strip_lab := paste0(cohort, "   (", lab, ")"), on = "cohort"]
pC <- ggplot(cf, aes(TotalExpr, RI_burden, colour = cohort)) +
  geom_point(size = .75, alpha = .55) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .55) +
  facet_wrap(~ strip_lab, scales = "free") +
  scale_colour_manual(values = COH_COL) +
  labs(x = "Global expression (mean log2 TPM across all genes)",
       y = "IR burden",
       title = "C  Global expression: the strongest confounder",
       subtitle = "stronger for low-expressed genes than high-expressed: partly technical") +
  TH + theme(legend.position = "none")

# ---------- D：LGG 免疫冷的直观证据 ----------
d <- sm[cohort == "LGG", .(sample, RI_burden, CD8_T, IFN_gamma_Ayers, Cytotoxic)]
d <- melt(d, id.vars = c("sample", "RI_burden"),
          variable.name = "panel", value.name = "score")
d[, q := cut(RI_burden, breaks = quantile(RI_burden, c(0, .25, .5, .75, 1)),
             labels = c("Q1", "Q2", "Q3", "Q4"),
             include.lowest = TRUE)]
setorderv(d, "RI_burden")
pD <- ggplot(d, aes(q, score)) +
  geom_boxplot(outlier.size = .3, linewidth = .28, width = .62,
               fill = "grey95", colour = "grey45") +
  geom_jitter(width = .12, size = .25, alpha = .25, colour = "#C0392B") +
  facet_wrap(~ panel, nrow = 1,
             labeller = as_labeller(c(CD8_T = "CD8 T",
                                      IFN_gamma_Ayers = "IFN-gamma",
                                      Cytotoxic = "Cytotoxic"))) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               colour = "#C0392B", linewidth = .45) +
  labs(x = "IR burden quartile", y = "Panel score (z, within cohort)",
       title = "D  LGG: colder with higher burden",
       subtitle = "red line joins group means") +
  TH

f <- (pA | pB) / (pC | pD) + theme(plot.margin = margin(4, 6, 4, 4))
ggsave(file.path(FIG, "Fig10_tme_confounded_coupling.png"), f,
       width = 200, height = 168, units = "mm", dpi = 600)
ggsave(file.path(FIG, "Fig10_tme_confounded_coupling.pdf"), f,
       width = 200, height = 168, units = "mm")
cat("Fig10 written\n")

# ---------- FigS2：GBM 同款箱线图（阴性对照）----------
g <- sm[cohort == "GBM", .(sample, RI_burden, CD8_T, IFN_gamma_Ayers, Cytotoxic)]
g <- melt(g, id.vars = c("sample", "RI_burden"),
          variable.name = "panel", value.name = "score")
g[, q := cut(RI_burden, breaks = quantile(RI_burden, c(0, .25, .5, .75, 1)),
             labels = c("Q1", "Q2", "Q3", "Q4"),
             include.lowest = TRUE)]
pS <- ggplot(g, aes(q, score)) +
  geom_boxplot(outlier.size = .3, linewidth = .28, width = .62,
               fill = "grey95", colour = "grey45") +
  geom_jitter(width = .12, size = .25, alpha = .25, colour = "#185FA5") +
  facet_wrap(~ panel, nrow = 1,
             labeller = as_labeller(c(CD8_T = "CD8 T",
                                      IFN_gamma_Ayers = "IFN-gamma",
                                      Cytotoxic = "Cytotoxic"))) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               colour = "#185FA5", linewidth = .45) +
  labs(x = "IR burden quartile", y = "Panel score (z, within cohort)",
       title = "GBM (n=150): no monotone immune gradient",
       subtitle = "unadjusted comparisons are null; a CD8 T/cytotoxic signal appears only after adjusting for global expression") +
  TH
ggsave(file.path(FIG, "FigS2_GBM_immune_null.png"), pS,
       width = 150, height = 62, units = "mm", dpi = 600)
ggsave(file.path(FIG, "FigS2_GBM_immune_null.pdf"), pS,
       width = 150, height = 62, units = "mm")
cat("FigS2 written\n")

# ---------- 供手稿摘数字 ----------
cat("\n=== LGG M7 (最保守模型) ===\n")
print(pc[cohort == "LGG" & model == M7 & panel %in% TME][order(rho),
      .(panel, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
        ci_hi = round(ci_hi, 3), p = signif(p, 3))])
cat("\n=== GBM M7 ===\n")
print(pc[cohort == "GBM" & model == M7 & panel %in% TME][order(rho),
      .(panel, rho = round(rho, 3), ci_lo = round(ci_lo, 3),
        ci_hi = round(ci_hi, 3), p = signif(p, 3))])
cat("\n=== LGG 四分位趋势检验（Spearman, burden vs score）===\n")
for (v in c("CD8_T", "IFN_gamma_Ayers", "Cytotoxic", "Endothelial")) {
  for (co in c("LGG", "GBM")) {
    x <- sm[cohort == co, ..v][[1]]; y <- sm[cohort == co, RI_burden]
    tt <- suppressWarnings(cor.test(x, y, method = "spearman"))
    cat(sprintf("%s %-16s rho=%+.3f p=%.2e n=%d\n", co, v,
                tt$estimate, tt$p.value, length(y)))
  }
}
