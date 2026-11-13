#!/usr/bin/env Rscript
# 43 — Fig9：泛癌剪接体耦合定位（把 GBM/LGG 放进 31 癌种坐标系）
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
OUT <- "output_pancan"; FIG <- "figures"
d <- fread(file.path(OUT, "pancan_sf_coupling.tsv"))
d[, cohort := sub("TCGA-", "", projectID)]
GLIO <- c("GBM", "LGG")

TH <- theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = .2),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        axis.text = element_text(size = 7),
        legend.position = "none")
pfmt <- function(p) ifelse(is.na(p), "NA",
                   ifelse(p < 1e-4, sprintf("p=%.0e", p), sprintf("p=%.3g", p)))
COLC <- c(`TRUE` = "#C0392B", `FALSE` = "grey70")

# 排序森林图
panel_forest <- function(pn, ttl, sub) {
  s <- d[panel == pn][order(rho)]
  s[, is_glio := cohort %in% GLIO]
  s[, cohort := factor(cohort, levels = cohort)]
  ggplot(s, aes(rho, cohort, colour = is_glio)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .18,
                   linewidth = .28, colour = "grey65") +
    geom_point(size = 1.7) +
    geom_text(data = s[(is_glio)], aes(x = ci_hi + .035,
              label = sprintf("%+.2f", rho)),
              size = 2.4, colour = "#C0392B", hjust = 0, fontface = "bold") +
    scale_colour_manual(values = COLC) +
    scale_x_continuous(expand = expansion(mult = c(.06, .18))) +
    labs(x = "Spearman rho with IR burden", y = NULL,
         title = ttl, subtitle = sub) + TH
}

pA <- panel_forest("SNRPC", "A  SNRPC: negative in 31/31 tumour types",
                   "median rho = -0.203; GBM -0.291, LGG -0.121")
pB <- panel_forest("Sm_core", "B  Sm core: GBM is the 3rd most negative of 31",
                   "GBM -0.305 vs LGG +0.105: the lineage reversal")
pC <- panel_forest("SR_proteins", "C  SR proteins: positive in 28/31",
                   "GBM +0.271, LGG +0.327")

# D：GBM/LGG 在泛癌分布中的位置（三个剪接体面板 + 两个对照）
key <- c("Sm_core", "SNRPC", "SR_proteins", "Proliferation",
         "Tumor_content", "NineteenComplex")
lab <- c(Sm_core = "Sm core", SNRPC = "SNRPC", SR_proteins = "SR proteins",
         Proliferation = "Proliferation", Tumor_content = "Tumour content",
         NineteenComplex = "NTC / Nineteen complex")
s <- d[panel %in% key]
s[, panel := factor(panel, levels = key)]
s[, labl := lab[as.character(panel)]]
gl <- s[cohort %in% GLIO]
bx <- s[, .(lo = min(rho), hi = max(rho), med = median(rho)), by = labl]
gl[, labl := lab[as.character(panel)]]
gl[, cohort := factor(cohort, levels = c("GBM", "LGG"))]

pD <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_boxplot(data = s, aes(labl, rho), width = .55, outlier.size = .6,
               fill = "grey92", colour = "grey55", linewidth = .25) +
  geom_point(data = gl, aes(labl, rho, colour = cohort, shape = cohort),
             size = 3.1, position = position_dodge(width = .38)) +
  geom_text(data = gl, aes(labl, rho, label = sprintf("%+.2f", rho),
            colour = cohort, group = cohort), vjust = -1.05, size = 2.4,
            position = position_dodge(width = .38), fontface = "bold") +
  scale_colour_manual(values = c(GBM = "#C0392B", LGG = "#2471A3")) +
  labs(x = NULL, y = "Spearman rho with IR burden",
       title = "D  Glioma positions within the pan-cancer distribution",
       subtitle = "boxes = 31 tumour types; GBM/LGG marked. Confounder panels stay at zero") +
  theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = .2),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = c(.86, .17),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3.2, "mm"),
        legend.background = element_rect(fill = "white", colour = "grey80"))

f <- (pA | pB) / (pC | pD)
ggsave(file.path(FIG, "Fig9_pancan_spliceosome_coupling.png"), f,
       width = 11.2, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(FIG, "Fig9_pancan_spliceosome_coupling.pdf"), f,
       width = 11.2, height = 8.2, bg = "white")
cat("[saved] Fig9_pancan_spliceosome_coupling\n")
