#!/usr/bin/env Rscript
# Step 23 — 手稿主图重制（Kahles 2018 双谱系口径）
#   Fig1 图谱 / Fig2 序列特征 / Fig3 连续轴 / Fig4 解耦+剪接体 / Fig7 非预后
#   （Fig5 机制锚定、Fig6 单细胞见 24_split_fig5_fig6.R，数据源独立未变）
suppressMessages({
  library(ggplot2); library(data.table); library(patchwork)
  library(survival); library(survminer); library(reshape2)
  library(dplyr); library(jsonlite); library(scales)
})

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
KM  <- paste0(OUT, "/output_km")
FIG <- paste0(OUT, "/figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

TH <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.2, colour = "grey35"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8))
CO <- c(GBM = "#C44E52", LGG = "#4C72B0")
COH <- c("GBM", "LGG")

save2 <- function(p, name, w, h) {
  ggsave(paste0(FIG, "/", name, ".png"), p, width = w, height = h, dpi = 300)
  pdf(paste0(FIG, "/", name, ".pdf"), width = w, height = h, onefile = TRUE)
  print(p); dev.off()
  cat(name, "done\n")
}

pfmt <- function(p) {                      # 向量化：可用于 geom_text 的 label
  ifelse(is.na(p), "NA",
  ifelse(p < 1e-4, sprintf("p = %.1e", p), sprintf("p = %.3g", p)))
}

## ===================== 共用数据 =====================
ev1  <- fread(paste0(KM, "/fig1_event_stats.tsv"))
sum1 <- fread(paste0(KM, "/fig1_summary.tsv"))
pcv  <- fread(paste0(KM, "/fig3_pcevr.tsv"))
sco  <- fread(paste0(KM, "/fig3_scores.tsv"))
ex4  <- fread(paste0(KM, "/fig4_expr.tsv"))
cdg  <- fromJSON(paste0(KM, "/cluster_diagnostics.json"))
g3   <- fread(paste0(KM, "/G3_km_intron_features.tsv"))

cf <- rbind(fread(paste0(KM, "/GBM/confound.tsv")), fread(paste0(KM, "/LGG/confound.tsv")))
sf <- rbind(fread(paste0(KM, "/GBM/sf_corr.tsv")), fread(paste0(KM, "/LGG/sf_corr.tsv")))
sv <- rbind(fread(paste0(KM, "/GBM/survival_merged.tsv")),
            fread(paste0(KM, "/LGG/survival_merged.tsv")))
# 队列归属：用 PCA 得分表里的 GBM 样本名（前 15 位）匹配
gbm_pat <- sco[cancer == "GBM", substr(sample, 1, 15)]
sv$cancer <- ifelse(substr(sv$sample, 1, 15) %in% gbm_pat, "GBM", "LGG")
sv$cancer <- factor(sv$cancer, levels = COH)
stopifnot(all(table(sv$cancer) > 40))

## ===================== Fig 1: 图谱 =====================
comp <- rbind(
  data.frame(cancer = factor(rep(COH, each = 3), levels = COH),
             stage  = factor(rep(c("IR events in atlas", "Detected in >=90% of tumours",
                                   "High-variability (SD >= 0.10)"), 2),
                             levels = c("IR events in atlas",
                                        "Detected in >=90% of tumours",
                                        "High-variability (SD >= 0.10)")),
             n = c(123390, 36682, 1522, 123390, 42990, 2216)))
lab1 <- data.frame(cancer = factor(rep(COH, each = 3), levels = COH),
                   stage = comp$stage[1:6 %% 3 == 1 | TRUE],
                   n = comp$n, y = comp$n * 1.25)
lab1$stage <- comp$stage
p1a <- ggplot(comp, aes(stage, n, fill = cancer)) +
  geom_col(position = position_dodge(.75), width = .68, colour = "grey25", linewidth = .2) +
  geom_text(aes(y = n * 1.3, label = format(n, big.mark = ",")),
            position = position_dodge(.75), size = 2.7, vjust = 0) +
  scale_fill_manual(values = CO, name = NULL) +
  scale_y_log10(expand = expansion(mult = c(0, .22)),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  labs(x = NULL, y = "IR events (log scale)",
       title = "A  IR event composition",
       subtitle = "SplAdder pan-cancer atlas; GBM n = 153, LGG n = 503") + TH +
  theme(axis.text.x = element_text(size = 8),
        legend.position = c(.86, .88))

p1b <- ggplot(sco, aes(RI_burden, fill = cancer)) +
  geom_histogram(aes(y = after_stat(density)), bins = 26,
                 colour = "white", linewidth = .2, alpha = .85) +
  geom_vline(aes(xintercept = mean(RI_burden)), linetype = 2, colour = "grey30") +
  scale_fill_manual(values = CO, name = NULL) +
  facet_wrap(~ cancer, scales = "free_y", nrow = 2) +
  labs(x = "Global IR burden (mean PSI of IR events detected in >=90% of tumours)",
       y = "Density",
       title = "B  Distribution of global IR burden",
       subtitle = sprintf("GBM mean %.3f (range %.3f-%.3f); LGG mean %.3f (range %.3f-%.3f)",
                          sum1[cancer == "GBM", burden_mean], sum1[cancer == "GBM", burden_min],
                          sum1[cancer == "GBM", burden_max], sum1[cancer == "LGG", burden_mean],
                          sum1[cancer == "LGG", burden_min], sum1[cancer == "LGG", burden_max])) +
  TH + theme(legend.position = "none")

set.seed(1)
evd <- ev1[!is.na(sd), .SD[sample(.N, min(.N, 30000))], by = cancer]
evd$hivar <- factor(ifelse(evd$hivar, "SD >= 0.10 (high-variability)", "SD < 0.10"),
                    levels = c("SD < 0.10", "SD >= 0.10 (high-variability)"))
p1c <- ggplot(evd, aes(detect, sd, colour = hivar)) +
  geom_point(size = .35, alpha = .35) +
  geom_hline(yintercept = 0.10, linetype = 3, colour = "grey40") +
  geom_vline(xintercept = 0.90, linetype = 3, colour = "grey40") +
  scale_colour_manual(values = c("SD < 0.10" = "grey70",
                                 "SD >= 0.10 (high-variability)" = "#E08214"),
                      name = NULL) +
  facet_wrap(~ cancer, nrow = 2) +
  labs(x = "Detection rate (fraction of tumours with non-missing PSI)",
       y = "SD of PSI across tumours",
       title = "C  Detection rate versus variability",
       subtitle = "random 30,000 events per cohort shown (events with <2 detections excluded); dashed lines = thresholds") +
  TH + theme(legend.position = c(.68, .85),
             legend.background = element_rect(fill = "white", colour = NA))

f1 <- p1a | (p1b / p1c) + plot_layout(widths = c(1.05, 1))
save2(f1, "Fig1_RI_atlas", 12.5, 8)

## ===================== Fig 2: 序列特征 =====================
g3$grp <- factor(g3$grp, levels = c("Spliced", "Variable", "Retained"))
kwp <- function(v, g, a = "Retained", b = "Spliced") {
  x <- v[g == a]; y <- v[g == b]
  suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value
}
boxp <- function(var, ttl, ylab, logy = FALSE, ylim = NULL) {
  d <- g3[, .(grp, v = get(var))]
  d <- d[!is.na(v)]
  p <- kwp(d$v, d$grp)
  gg <- ggplot(d, aes(grp, v, fill = grp)) +
    geom_boxplot(outlier.shape = NA, width = .62, linewidth = .3) +
    scale_fill_manual(values = c(Spliced = "#9FB6CD", Variable = "#EDC775",
                                 Retained = "#C44E52"), guide = "none") +
    labs(x = NULL, y = ylab, title = ttl,
         subtitle = sprintf("median %.3g / %.3g / %.3g;  Wilcoxon %s",
                            median(d$v[d$grp == "Spliced"]),
                            median(d$v[d$grp == "Variable"]),
                            median(d$v[d$grp == "Retained"]), pfmt(p))) + TH
  if (logy) gg <- gg + scale_y_log10()
  if (!is.null(ylim)) gg <- gg + coord_cartesian(ylim = ylim)
  gg + theme(axis.text.x = element_text(size = 8.5))
}
p2a <- boxp("ilen", "A  Intron length", "Intron length (bp)", logy = TRUE)
p2b <- boxp("gc", "B  GC content", "GC content", ylim = c(0.3, 0.8))
p2c <- boxp("y20", "C  Polypyrimidine-tract strength", "y20 (poly-Y score)", ylim = c(0.35, 0.9))
p2d <- boxp("e2len", "D  Downstream exon length", "Downstream exon length (bp)", logy = TRUE)

# E: 长度匹配后（19 个 ilen 分位层内等数抽样）
d2 <- g3[grp %in% c("Retained", "Spliced")]
nb <- 19
d2$bin <- as.integer(cut(d2$ilen,
                         breaks = quantile(d2$ilen, probs = seq(0, 1, length.out = nb + 1)),
                         include.lowest = TRUE))
set.seed(7)
mt <- rbindlist(lapply(split(d2, d2$bin), function(s) {
  nr <- sum(s$grp == "Retained"); ns <- sum(s$grp == "Spliced")
  if (nr < 2 || ns < 2) return(NULL)
  n <- min(nr, ns)
  rbind(s[grp == "Retained"][sample(nr, n)], s[grp == "Spliced"][sample(ns, n)])
}))
mtl <- melt(mt[, .(grp, gc, y20)], id.vars = "grp")
mtl$variable <- factor(mtl$variable, levels = c("gc", "y20"),
                       labels = c("GC content", "poly-Y (y20)"))
ps <- sapply(c("gc", "y20"), function(v)
  suppressWarnings(wilcox.test(mt[grp == "Retained", get(v)],
                               mt[grp == "Spliced", get(v)], exact = FALSE))$p.value)
mtl$grp <- factor(mtl$grp, levels = c("Spliced", "Retained"))
p2e <- ggplot(mtl, aes(grp, value, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = .55, linewidth = .3) +
  scale_fill_manual(values = c(Spliced = "#9FB6CD", Retained = "#C44E52"), guide = "none") +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = NULL, y = "Value", title = "E  Length-matched comparison",
       subtitle = sprintf("n = %d per group after matching on intron-length ventiles; GC %s, y20 %s",
                          unique(table(mt$grp)), pfmt(ps[1]), pfmt(ps[2]))) + TH

f2 <- (p2a | p2b | p2c) / (p2d | p2e | plot_spacer())
save2(f2, "Fig2_sequence_features", 12.5, 7.5)

## ===================== Fig 3: 连续轴 =====================
pcv$PC <- factor(pcv$PC, levels = paste0("PC", 1:10))
pcv$cancer <- factor(pcv$cancer, levels = COH)
p3a <- ggplot(pcv, aes(PC, var, fill = cancer)) +
  geom_col(position = position_dodge(.75), width = .7, colour = "grey25", linewidth = .15) +
  scale_fill_manual(values = CO, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Variance explained",
       title = "A  Principal components of high-variability IR events",
       subtitle = sprintf("GBM %d events (n = 153); LGG %d events (n = 503)",
                          sum1[cancer == "GBM", n_hivar], sum1[cancer == "LGG", n_hivar])) + TH +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = c(.85, .78))

rho_lab <- sco[, .(rho = cor(PC1, RI_burden, method = "spearman")), by = cancer]
setkey(rho_lab, cancer); setkey(sco, cancer)
p3b <- ggplot(sco, aes(PC1, RI_burden, colour = cancer)) +
  geom_point(size = 1.3, alpha = .8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  scale_colour_manual(values = CO, name = NULL) +
  facet_wrap(~ cancer, scales = "free", nrow = 2) +
  geom_text(data = rho_lab, aes(x = Inf, y = -Inf,
                                label = sprintf("Spearman rho = %.3f", rho)),
            hjust = 1.1, vjust = -1.2, size = 3.2, colour = "grey25", inherit.aes = FALSE) +
  labs(x = "PC1 of high-variability IR events", y = "Global IR burden",
       title = "B  PC1 is essentially the IR-burden axis",
       subtitle = "single dominant axis; no separated islands") + TH +
  theme(legend.position = "none")

pac_dt <- rbindlist(lapply(COH, function(c) {
  x <- cdg[[c]]$PAC
  data.frame(cancer = c, k = as.integer(names(x)), pac = unlist(x))
}))
sil_dt <- data.frame(cancer = COH,
                     sil = sapply(COH, function(c) cdg[[c]]$silhouette_k2))
dau_dt <- rbindlist(lapply(COH, function(c) {
  x <- cdg[[c]]$delta_AUC
  x <- x[!sapply(x, is.null)]
  data.frame(cancer = c, k = as.integer(names(x)), dauc = unlist(x))
}))
bic_dt <- rbindlist(lapply(COH, function(c) {
  x <- cdg[[c]]$BIC
  data.frame(cancer = c, k = factor(paste0("k=", names(x)), levels = c("k=1", "k=2")),
             bic = unlist(x))
}))
pac_dt$cancer <- factor(pac_dt$cancer, levels = COH)
dau_dt$cancer <- factor(dau_dt$cancer, levels = COH)
sil_dt$cancer <- factor(sil_dt$cancer, levels = COH)
bic_dt$cancer <- factor(bic_dt$cancer, levels = COH)

p3c1 <- ggplot(pac_dt, aes(k, pac, colour = cancer, group = cancer)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = CO, name = NULL) +
  scale_x_continuous(breaks = 2:6) +
  labs(x = "k", y = "PAC (corrected)", title = "C  Clustering diagnostics",
       subtitle = "low PAC = reproducible split, not separation; read with silhouette") + TH +
  theme(legend.position = c(.8, .82))
p3c2 <- ggplot(sil_dt, aes(cancer, sil, fill = cancer)) +
  geom_col(width = .55, colour = "grey25", linewidth = .2) +
  geom_text(aes(y = sil + .012, label = sprintf("%.3f", sil)), size = 2.9) +
  scale_fill_manual(values = CO, guide = "none") +
  coord_cartesian(ylim = c(0, .2)) +
  labs(x = NULL, y = "Silhouette (k = 2)",
       title = "D  Low silhouette at k = 2",
       subtitle = "well below the ~0.5 expected for real clusters") + TH
p3c3 <- ggplot(dau_dt, aes(k, dauc, colour = cancer, group = cancer)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  scale_colour_manual(values = CO, name = NULL) +
  scale_x_continuous(breaks = 3:6) +
  labs(x = "k", y = expression(Delta * "AUC"),
       title = "E  No elbow in the consensus CDF",
       subtitle = "gains beyond k = 2 are marginal") + TH +
  theme(legend.position = "none")
p3c4 <- ggplot(bic_dt, aes(cancer, bic, fill = k)) +
  geom_col(position = position_dodge(.7), width = .62, colour = "grey25", linewidth = .2) +
  scale_fill_manual(values = c("k=1" = "#8DA0CB", "k=2" = "#FC8D62"), name = NULL) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(x = NULL, y = "GMM BIC (log scale)",
       title = "F  One-component model preferred",
       subtitle = "lower BIC = better; k = 1 wins in both cohorts") + TH +
  theme(legend.position = c(.85, .88))

f3 <- (p3a | p3b) / (p3c1 | p3c2 | p3c3 | p3c4)
save2(f3, "Fig3_continuous_axis", 13, 9)

## ===================== Fig 4: 解耦 + 剪接体 =====================
cf$cancer <- factor(rep(COH, each = nrow(cf) / 2), levels = COH)
cfc <- cf[`gene/panel` %in% c("Proliferation", "TumorPurity_marker", "Neuronal_normal")]
cfc$panel <- factor(cfc$`gene/panel`,
                    levels = c("Proliferation", "TumorPurity_marker", "Neuronal_normal"),
                    labels = c("Proliferation", "Tumour-content markers", "Neuronal / normal tissue"))
p4a <- ggplot(cfc, aes(panel, rho, fill = cancer)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_col(position = position_dodge(.72), width = .62, colour = "grey25", linewidth = .2) +
  geom_text(aes(y = rho + sign(rho) * .03, label = sprintf("%.3f\n%s", rho, pfmt(p))),
            position = position_dodge(.72), size = 2.4, lineheight = .85) +
  scale_fill_manual(values = CO, name = NULL) +
  coord_cartesian(ylim = c(-.28, .28)) +
  labs(x = NULL, y = "Spearman rho with IR burden",
       title = "A  IR burden is decoupled from proliferation and tumour content",
       subtitle = "GBM n = 151, LGG n = 503") + TH +
  theme(axis.text.x = element_text(size = 8), legend.position = c(.85, .88))

sfx <- cf[`gene/panel` %in% c("Sm_core", "SR_proteins", "NineTeen_Bact")]
sfx$panel <- factor(sfx$`gene/panel`, levels = c("Sm_core", "SR_proteins", "NineTeen_Bact"),
                    labels = c("Sm core (assembly core)", "SR proteins", "Nineteen / Bact complex"))
p4b <- ggplot(sfx, aes(panel, rho, fill = cancer)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_col(position = position_dodge(.72), width = .62, colour = "grey25", linewidth = .2) +
  geom_text(aes(y = rho + sign(rho) * .04, label = sprintf("%.3f\n%s", rho, pfmt(p))),
            position = position_dodge(.72), size = 2.4, lineheight = .85) +
  scale_fill_manual(values = CO, name = NULL) +
  coord_cartesian(ylim = c(-.55, .8)) +
  labs(x = NULL, y = "Spearman rho (sub-complex mean expression)",
       title = "B  Spliceosome sub-complexes: Sm-core association reverses in LGG",
       subtitle = "sub-complex averaging avoids per-gene pseudo-replication") + TH +
  theme(axis.text.x = element_text(size = 8), legend.position = c(.85, .88))

ex4$cancer <- factor(ex4$cancer, levels = COH)
rl <- ex4[, .(rho = cor(SNRPC_z, RI_burden, method = "spearman"),
              p = cor.test(SNRPC_z, RI_burden, method = "spearman")$p.value), by = cancer]
p4c <- ggplot(ex4, aes(SNRPC_z, RI_burden, colour = cancer)) +
  geom_point(size = 1.1, alpha = .75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  scale_colour_manual(values = CO, name = NULL) +
  facet_wrap(~ cancer, scales = "free", nrow = 2) +
  geom_text(data = rl, aes(x = -Inf, y = -Inf,
                           label = sprintf("rho = %.3f, %s", rho, pfmt(p))),
            hjust = -0.12, vjust = -1.4, size = 3.2, colour = "grey25", inherit.aes = FALSE) +
  labs(x = "SNRPC expression (z-score across tumours)", y = "Global IR burden",
       title = "C  SNRPC versus IR burden",
       subtitle = "each point = one tumour") + TH + theme(legend.position = "none")

sf$gene <- factor(sf$gene, levels = c("SNRPC", "SNRPB", "SNRPD1", "SNRPD2", "SNRPD3",
                                      "SNRPE", "SNRPF", "SNRPG", "SRSF1", "SRSF2", "SRSF3",
                                      "SRSF6", "SRSF7", "SRSF8"))
sf$class <- ifelse(as.character(sf$gene) %in%
                     c("SNRPC", "SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"),
                   "Sm core / SNRPC", "SR proteins")
sf$cancer <- factor(rep(COH, each = nrow(sf) / 2), levels = COH)
p4d <- ggplot(sf, aes(gene, rho, fill = cancer)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_col(position = position_dodge(.75), width = .66, colour = "grey25", linewidth = .15) +
  facet_wrap(~ class, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = CO, name = NULL) +
  geom_text(data = sf[cancer == "GBM"], aes(y = rho - .045,
            label = ifelse(p < .001, "***", ifelse(p < .01, "**",
                    ifelse(p < .05, "*", "ns")))), size = 2.4, colour = "grey35") +
  labs(x = NULL, y = "Spearman rho with IR burden",
       title = "D  Per-gene correlations",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001 (GBM significance marked on bars)") + TH +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top")

f4 <- (p4a | p4b) / (p4c | p4d)
save2(f4, "Fig4_decoupling_spliceosome", 12.5, 8.5)

## ===================== Fig 7: 非预后 =====================
sv$burden_group <- factor(sv$burden_group, levels = c("RI-low", "RI-high"))
kmplot <- function(c) {
  d <- sv[cancer == c & !is.na(OS.time) & !is.na(OS) & OS.time > 0]
  fit <- survfit(Surv(OS.time, OS) ~ burden_group, data = d)
  lr <- survdiff(Surv(OS.time, OS) ~ burden_group, data = d)
  p <- 1 - pchisq(lr$chisq, 1)
  g <- ggsurvplot(fit, data = d, pval = FALSE, risk.table = FALSE,
                  legend.title = "IR burden", legend.labs = c("low", "high"),
                  palette = c("#4C72B0", "#C44E52"), xlab = "Time (days)",
                  ylab = "Overall survival", surv.median.line = "none",
                  ggtheme = TH, font.x = 9, font.y = 9, font.tickslab = 8,
                  title = sprintf("%s (n = %d)", c, nrow(d)),
                  font.title = c(11, "bold", "black"),
                  break.time.by = ifelse(c == "GBM", 500, 1000))
  g$plot <- g$plot + annotate("text", x = Inf, y = Inf,
                              label = sprintf("log-rank %s", pfmt(p)),
                              hjust = 1.15, vjust = 1.8, size = 3.2, colour = "grey30")
  g$plot
}
p7a1 <- kmplot("GBM"); p7a2 <- kmplot("LGG")
p7a <- (p7a1 | p7a2) + plot_annotation(
  title = "A  IR burden does not stratify overall survival",
  subtitle = "median split of IR burden; TCGA-GBM and TCGA-LGG",
  theme = theme(plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 8.2, colour = "grey35")))

# B: 森林图（burden 每 SD vs 年龄每 10 岁，阳性对照）
fit_rows <- list()
for (c in COH) {
  d <- sv[cancer == c & !is.na(OS.time) & !is.na(OS) & OS.time > 0]
  bz <- as.numeric(scale(d$RI_burden))
  az <- d$age / 10          # 每 10 岁，不做标准化
  for (ep in c("OS", "DSS", "PFI")) {
    tt <- d[[paste0(ep, ".time")]]; ee <- d[[ep]]
    ok <- !is.na(tt) & !is.na(ee) & tt > 0
    m <- coxph(Surv(tt[ok], ee[ok]) ~ bz[ok])
    s <- summary(m)
    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      cohort = c, term = sprintf("IR burden / SD (%s)", ep),
      HR = s$conf.int[1], lo = s$conf.int[3], hi = s$conf.int[4],
      p = s$coefficients[5], ctrl = FALSE)
  }
  m <- coxph(Surv(d$OS.time, d$OS) ~ az)
  s <- summary(m)
  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    cohort = c, term = "Age / 10 years (OS)", HR = s$conf.int[1],
    lo = s$conf.int[3], hi = s$conf.int[4], p = s$coefficients[5], ctrl = TRUE)
}
fr <- rbindlist(fit_rows)
fr$term <- factor(fr$term, levels = rev(c("IR burden / SD (OS)", "IR burden / SD (DSS)",
                                          "IR burden / SD (PFI)", "Age / 10 years (OS)")))
fr$cohort <- factor(fr$cohort, levels = COH)
p7b <- ggplot(fr, aes(term, HR, colour = cohort, shape = ctrl)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .12,
                position = position_dodge(.5), linewidth = .4) +
  geom_point(position = position_dodge(.5), size = 2.2, fill = "white") +
  scale_colour_manual(values = CO, name = NULL) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("IR burden", "Age (positive control)"), name = NULL) +
  scale_x_discrete(labels = function(x) sub("IR burden / SD", "IR burden/SD", x)) +
  coord_flip(ylim = c(0.55, 3.2)) +
  geom_text(aes(y = hi + .06, label = sprintf("%.2f (%s)", HR, pfmt(p))),
            position = position_dodge(.5), size = 2.5, hjust = 0) +
  labs(x = NULL, y = "Hazard ratio (95% CI)",
       title = "B  Age positive control: the null is not a power artefact",
       subtitle = "filled circles = IR burden; triangles = age-at-diagnosis (per 10 years)") + TH +
  theme(legend.position = "bottom", legend.box = "horizontal",
        axis.text.y = element_text(size = 8))

cg <- fread(paste0(OUT, "/output/cgga_cox_univariate.tsv"))
cg$subset_en <- c("全队列" = "Whole cohort (n = 970)",
                  "WHO IV (GBM)" = "WHO IV (GBM)",
                  "WHO III" = "WHO III", "WHO II" = "WHO II",
                  "IDH-wildtype" = "IDH-wildtype", "IDH-mutant" = "IDH-mutant",
                  "IDH-wt WHO IV" = "IDH-wt WHO IV",
                  "1p19q codel" = "1p/19q co-deleted",
                  "1p19q non-codel" = "1p/19q non-co-deleted")[cg$subset]
cg$subset_en <- factor(cg$subset_en, levels = rev(c(
  "Whole cohort (n = 970)", "WHO IV (GBM)", "IDH-wt WHO IV", "IDH-wildtype",
  "WHO III", "WHO II", "IDH-mutant", "1p/19q co-deleted", "1p/19q non-co-deleted")))
p7c <- ggplot(cg, aes(subset_en, HR_per_SD)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12, linewidth = .35) +
  geom_point(size = 2.4, colour = "#7B6888", fill = "white", shape = 21) +
  coord_flip(ylim = c(0.55, 1.6)) +
  geom_text(aes(y = upper + .03, label = sprintf("%.2f (%s)", HR_per_SD, pfmt(p))),
            size = 2.5, hjust = 0) +
  labs(x = NULL, y = "Hazard ratio per +1 SD (95% CI)",
       title = "C  CGGA: association confined to grade / IDH strata",
       subtitle = "expression-level RI-axis surrogate, CGGA n = 970; adjusted HR 0.97 (p = 0.54)") + TH +
  theme(axis.text.y = element_text(size = 8))

f7 <- (p7a) / (p7b | p7c) + plot_layout(heights = c(1, 1.15))
save2(f7, "Fig7_prognostic_null", 13, 11)

cat("\nFig1/2/3/4/7 完成 -> ", FIG, "\n")
