#!/usr/bin/env Rscript
# Step 24 — 拆分并重存 Fig5（机制锚定，GSE227187）与 Fig6（单细胞定位，GSE182109）
# 数据源独立于 Kahles 图谱，数字未变；仅按手稿新编号拆成两张图。
suppressMessages({
  library(ggplot2); library(patchwork); library(reshape2)
})

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
FIG <- paste0(OUT, "/figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

TH <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.2, colour = "grey35"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8))

save2 <- function(p, name, w, h) {
  ggsave(paste0(FIG, "/", name, ".png"), p, width = w, height = h, dpi = 300)
  pdf(paste0(FIG, "/", name, ".pdf"), width = w, height = h, onefile = TRUE)
  print(p); dev.off()
  cat(name, "done\n")
}

## ---------- Fig 5: SNRPC 敲低 → IR + OXPHOS 抑制（GSE227187，乳腺癌细胞系） ----------
ev <- data.frame(type = c("SE", "A5SS", "A3SS", "MXE", "RI"),
                 up = c(2485, 802, 940, 1024, 794),
                 down = c(2764, 670, 776, 1176, 606))
evl <- melt(ev, id.vars = "type", variable.name = "dir", value.name = "n")
evl$dir <- factor(evl$dir, levels = c("down", "up"),
                  labels = c("inclusion down", "inclusion up"))
p5a <- ggplot(evl, aes(type, n, fill = dir)) + geom_col(position = "dodge", width = .72) +
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
p5b <- ggplot(ox, aes(med, set)) +
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_col(aes(fill = pct), width = .68) +
  geom_text(aes(label = sprintf("%d%% down, n=%d", pct, n)), hjust = -0.08, size = 2.9) +
  scale_fill_gradient(low = "#F2C1C1", high = "#A32D2D", name = "% down") +
  coord_cartesian(xlim = c(-0.32, 0.06)) +
  labs(x = "Median log2 fold change (SNRPC KD)", y = NULL,
       title = "B  OXPHOS genes are coordinately down-regulated",
       subtitle = "167 OXPHOS genes, 77% down, Wilcoxon p = 5.8e-17") + TH +
  theme(legend.position = "none")

p5c <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  annotate("rect", xmin = .5, xmax = 5.5, ymin = .55, ymax = .95,
           fill = "#EAF3DE", colour = "#639922", linewidth = .4) +
  annotate("rect", xmin = .5, xmax = 5.5, ymin = .1, ymax = .5,
           fill = "#F1EFE8", colour = "#B4B2A9", linewidth = .4) +
  annotate("text", x = 3, y = .82, size = 3.4, fontface = "bold", colour = "#27500A",
           label = "21 / 30 genes with increased RI show reduced mRNA (70%)") +
  annotate("text", x = 3, y = .66, size = 3,
           label = "host-gene mRNA loss concordant with retention gain") +
  annotate("text", x = 3, y = .30, size = 3.1, colour = "#444441",
           label = "Representative: UQCRC1 RI +0.22 / log2FC -0.37") +
  annotate("text", x = 3, y = .16, size = 3.1, colour = "#444441",
           label = "COX7C RI +0.33 / log2FC -0.15;  MRS2 RI +0.65") +
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 1)) + theme_void() +
  labs(title = "C  Intron retention is coupled to transcript loss",
       subtitle = "Breast-cancer cell line; interpreted as cross-cancer mechanistic corroboration") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.2, colour = "grey40"))

f5 <- (p5a | p5b) / p5c
save2(f5, "Fig5_mechanism_anchor", 12, 8)

## ---------- Fig 6: 单细胞定位（GSE182109） ----------
cta <- read.delim(paste0(OUT, "/output/sc_celltype_RIaxis.tsv"), row.names = 1)
cta$cell_type <- rownames(cta)
cta$is_tumour <- cta$cell_type == "Tumor"

## A: pos-only（描述性口径，无内部对照）
ctaA <- cta
ctaA$cell_type <- factor(ctaA$cell_type, levels = ctaA$cell_type[order(ctaA$RIaxis_pos)])
p6a <- ggplot(ctaA, aes(RIaxis_pos, cell_type, fill = is_tumour)) +
  geom_col(width = .7) +
  scale_fill_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#B4B2A9"), guide = "none") +
  geom_text(aes(label = sprintf("n = %s", format(n_cells, big.mark = ","))),
            hjust = -0.05, size = 2.6) +
  coord_cartesian(xlim = c(0, max(ctaA$RIaxis_pos) * 1.45)) +
  labs(x = "Mean expression of RI-axis positive genes (log1p CP10K)", y = NULL,
       title = "A  Positive-genes-only metric (descriptive)",
       subtitle = "no internal control: tracks per-cell output, overstates the gap ~5x") + TH

## B: pos - neg（内部对照口径）
ctaB <- cta
ctaB$cell_type <- factor(ctaB$cell_type, levels = ctaB$cell_type[order(ctaB$RIaxis_score)])
p6b <- ggplot(ctaB, aes(RIaxis_score, cell_type, fill = is_tumour)) +
  geom_col(width = .7) +
  scale_fill_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#B4B2A9"), guide = "none") +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey45") +
  coord_cartesian(xlim = c(min(ctaB$RIaxis_score) * 1.25,
                           max(ctaB$RIaxis_score) * 1.35)) +
  labs(x = "RI-axis score (positive minus negative genes)", y = NULL,
       title = "B  Internal-control metric",
       subtitle = "with output controlled, compartment differences are small and not tumour-maximal") + TH

p6 <- (p6a | p6b) +
  plot_annotation(
    title = "The IR axis is anchored in malignant cells — but only modestly",
    subtitle = "GSE182109, 264,523 cells. Panel A has no internal control and must not be read as a compartment effect; the quantified tumour-immune difference (+0.026, p = 0.0058) comes from within-sample pairing, not from these marginal means",
    theme = theme(plot.title = element_text(size = 11.5),
                  plot.subtitle = element_text(size = 8.5, colour = "grey35")))
save2(p6, "Fig6_single_cell", 12.4, 5)

cat("\nFig5/Fig6 完成 -> ", FIG, "\n")
