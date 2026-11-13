#!/usr/bin/env Rscript
# G1 / G2 结果出图
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})
OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output"
FIG <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

TH <- theme_bw() + theme(plot.title = element_text(size = 12, face = "bold"),
                         axis.title = element_text(size = 11),
                         axis.text = element_text(size = 9),
                         strip.text = element_text(size = 10))

save_both <- function(g, name, w, h) {
  png(file.path(FIG, paste0(name, ".png")), width = w * 300, height = h * 300, res = 300)
  print(g); dev.off()
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h)
  print(g); dev.off()
}

## ============ Fig8: G1 negative control ============
ctl <- read.delim(file.path(OUT, "G1_negative_control_eventtypes.tsv"), check.names = FALSE)
ctl_main <- ctl[grepl("（全部）", ctl$type), ]
ctl_main$etype <- sub("（全部）", "", ctl_main$type)
ctl_main$etype <- factor(ctl_main$etype, levels = c("IR", "ES", "A3", "A5", "ME"))
ctl_main$sig <- ifelse(ctl_main$p_perm < 0.05, "perm p<0.05", "n.s.")

pA <- ggplot(ctl_main, aes(x = etype, y = frac_neg * 100, fill = sig)) +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = (null_mean - 1.96 * null_sd) * 100,
                     ymax = (null_mean + 1.96 * null_sd) * 100),
                 colour = "grey20", linewidth = 0.4) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.6, size = 2.8) +
  scale_fill_manual(values = c("perm p<0.05" = "#D85A30", "n.s." = "#B4B2A9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, 85)) +
  labs(x = "Event type", y = "% events with negative rho (PSI vs host-gene mRNA)",
       title = "A  Negative-coupling rate by event type",
       subtitle = "Dashed line = 50% chance; error bars = 95% permutation null",
       fill = "") + TH +
  theme(legend.position = c(0.75, 0.88), legend.background = element_blank())

rho <- read.delim(file.path(OUT, "G1_rho_distribution.tsv"))
rho$type <- factor(rho$type, levels = c("IR", "ES", "A5"))
labdat <- rho %>% group_by(type) %>%
  summarise(med = median(rho), n = n(), neg = mean(rho < 0) * 100, .groups = "drop")

pB <- ggplot(rho, aes(x = rho, colour = type, fill = type)) +
  geom_density(alpha = 0.25, linewidth = 0.6, adjust = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("IR" = "#D85A30", "ES" = "#378ADD", "A5" = "#888780")) +
  scale_fill_manual(values = c("IR" = "#D85A30", "ES" = "#378ADD", "A5" = "#888780")) +
  annotate("text", x = -0.55, y = c(4.2, 3.7, 3.2), size = 2.9, colour = "#333333",
           label = sprintf("%s: median rho %+.3f, %.1f%% negative (n=%s)",
                           labdat$type, labdat$med, labdat$neg, format(labdat$n, big.mark = ","))) +
  labs(x = expression(rho ~ "(intron-retention PSI vs host-gene mRNA)"),
       y = "Density", colour = "", fill = "",
       title = "B  Distribution of coupling coefficients",
       subtitle = "IR shifted furthest to the left; A5 centred at zero") + TH +
  theme(legend.position = "none")

fig8 <- ggarrange(pA, pB, ncol = 2, nrow = 1, widths = c(1, 1.15))
save_both(fig8, "Fig8_G1_negative_control", 11.5, 5.0)

## ============ Fig9: G2 spliceosome machinery ============
sf <- read.delim(file.path(OUT, "G2_sf_vs_burden.tsv"))
grp <- read.delim(file.path(OUT, "G2b_complex_groups.tsv"))
sc <- read.delim(file.path(OUT, "G2_scatter_data.tsv"))

# A: SNRPC scatter
sp <- cor.test(sc$SNRPC, sc$RI_burden, method = "spearman")
pA2 <- ggplot(sc, aes(x = SNRPC, y = RI_burden)) +
  geom_point(size = 2.0, alpha = 0.75, colour = "#185FA5") +
  geom_smooth(method = "lm", se = TRUE, colour = "#D85A30", fill = "#F5C4B3", linewidth = 0.7) +
  annotate("text", x = min(sc$SNRPC, na.rm = TRUE), y = max(sc$RI_burden, na.rm = TRUE),
           hjust = 0, vjust = 1, size = 3.2,
           label = sprintf("Spearman rho = %+.3f\np = %.2g  (n = %d)",
                           sp$estimate, sp$p.value, nrow(sc))) +
  labs(x = "SNRPC expression (log2)", y = "Global RI burden",
       title = "A  SNRPC (U1-C) vs intron-retention burden",
       subtitle = "TCGA-GBM, n = 151 primary tumours") + TH

# B: Sm core + NTC lollipop
sm <- c("SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG")
ntc <- c("ISY1", "ZMAT2", "CDC5L", "BCAS2", "BUD31", "CRNKL1", "AQR", "XAB2",
         "SYF2", "CWC15", "PLRG1", "PRPF19")
sel <- sf %>% filter(gene %in% c(sm, ntc)) %>%
  mutate(complex = ifelse(gene %in% sm, "Sm core (shared ring)", "NineTeen / Bact")) %>%
  arrange(rho) %>% mutate(gene = factor(gene, levels = gene))
pB2 <- ggplot(sel, aes(x = gene, y = rho, fill = complex)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c("Sm core (shared ring)" = "#534AB7",
                               "NineTeen / Bact" = "#1D9E75")) +
  coord_flip() + labs(x = "", y = "Spearman rho with RI burden", fill = "",
                      title = "B  Assembly / catalysis core: uniformly negative",
                      subtitle = "Sm core 7/7 negative, permutation p = 0.022") +
  TH + theme(legend.position = c(0.70, 0.16),
             legend.background = element_blank())

# C: complex-level forest
EN <- c("Sm 核心（U1/U2/U4/U5 共用）" = "Sm core (shared ring)",
        "U1 特异亚基" = "U1-specific subunits",
        "U2 / 3' 剪接位点识别" = "U2 / 3' splice site",
        "U4/U6.U5 tri-snRNP" = "U4/U6.U5 tri-snRNP",
        "NineTeen / Bact 复合体" = "NineTeen / Bact",
        "SR 蛋白家族" = "SR proteins",
        "hnRNP 家族" = "hnRNP family",
        "3' 端加工 / 输出" = "3'-end processing / export")
grp$lab <- ifelse(grp$group %in% names(EN), EN[grp$group], grp$group)
grp$dir <- ifelse(grp$median_rho < 0, "negative", "positive")
grp$lab <- factor(grp$lab, levels = grp$lab[order(grp$median_rho)])
pC2 <- ggplot(grp, aes(x = lab, y = median_rho)) +
  geom_segment(aes(xend = lab, y = 0, yend = median_rho, colour = dir),
               linewidth = 0.9) +
  geom_point(aes(colour = dir), size = 3.4) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%d/%d", n_neg, n),
                hjust = ifelse(median_rho < 0, 1.35, -0.35)),
            size = 2.7, colour = "grey30") +
  scale_colour_manual(values = c("negative" = "#185FA5", "positive" = "#D85A30"),
                      labels = c("negative" = "towards negative rho",
                                 "positive" = "towards positive rho")) +
  coord_flip(ylim = range(grp$median_rho) * 1.35) +
  labs(x = "", y = "Median rho with RI burden", colour = "",
       title = "C  Direction of association by sub-complex",
       subtitle = "Labels = members with negative rho / total") + TH +
  theme(legend.position = c(0.74, 0.13), legend.background = element_blank(),
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8))

fig9 <- ggarrange(pA2, ggarrange(pB2, pC2, ncol = 2, widths = c(1, 1.05)),
                  ncol = 1, heights = c(1, 1.05))
save_both(fig9, "Fig9_G2_spliceosome_core", 11.5, 9.2)

cat("图已输出：", FIG, "\n")
