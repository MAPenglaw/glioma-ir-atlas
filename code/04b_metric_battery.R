#!/usr/bin/env Rscript
# Step 4b — RI 指标电池：穷举式检验「是否存在任何有预后价值的 RI 衍生指标」
# 目的：不是挑选显著指标（那是 p-hacking），而是判定「RI 轴整体是否与预后相关」

suppressMessages({ library(survival); library(data.table) })

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
tab <- readRDS(paste0(OUT, "/output/survival_input.rds"))
## PC1-3 在 annotated 表中，补进来
ann <- read.delim(paste0(OUT, "/output/sample_RI_group_annotated.tsv"),
                  row.names = 1, check.names = FALSE)
pc_cols <- intersect(c("RI_PC1", "RI_PC2", "RI_PC3"), colnames(ann))
tab <- cbind(tab, ann[rownames(tab), pc_cols, drop = FALSE])

## 追加几种常见的 RI 衍生指标
b <- as.numeric(tab$RI_burden)
tab$burden_hivar  <- as.numeric(tab$RI_burden_hivar)
tab$burden_median <- as.numeric(tab$RI_burden_median)
## PSI>0.5 的 IR 事件占比（“滞留事件比例”口径）
psi50 <- readRDS(paste0(OUT, "/output/ri_burden_extra.rds"))
tab$frac_psi_gt_05 <- psi50[rownames(tab), "frac_psi_gt_05"]
tab$n_psi_gt_05    <- psi50[rownames(tab), "n_psi_gt_05"]
tab$burden_top200  <- psi50[rownames(tab), "burden_top200"]

VARS <- c("RI_burden", "burden_hivar", "burden_median", "frac_psi_gt_05",
          "n_psi_gt_05", "burden_top200", "RI_PC1", "RI_PC2", "RI_PC3",
          "U1_snRNP", "Spliceosome_core", "Proliferation",
          "Glioma_StemCell", "IFN_gamma_signature", "Treg", "CD8_T_cell")

out <- data.frame()
for (v in VARS) {
  if (!v %in% colnames(tab)) { cat("  [跳过] 缺列:", v, "\n"); next }
  x <- suppressWarnings(as.numeric(tab[[v]]))
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) { cat("  [跳过] 无变异:", v, "\n"); next }
  tab[[v]] <- x
  tab[[paste0("z_", v)]] <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  for (ep in c("OS", "PFI")) {
    d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]) &
             !is.na(tab[[paste0("z_", v)]]), ]
    if (nrow(d) < 60) next
    m <- coxph(as.formula(paste0("Surv(", ep, ".time, ", ep, ") ~ z_", v)), data = d)
    s <- summary(m)
    out <- rbind(out, data.frame(
      endpoint = ep, metric = v, n = m$n, events = m$nevent,
      HR_per_SD = s$coefficients[1, "exp(coef)"],
      lower = s$conf.int[1, "lower .95"], upper = s$conf.int[1, "upper .95"],
      p = s$coefficients[1, "Pr(>|z|)"]))
  }
}
out <- out[order(out$endpoint, out$p), ]
out$p_adj_within_endpoint <- ave(out$p, out$endpoint, FUN = function(x) p.adjust(x, "BH"))
cat("\n================ 指标电池：每 +1 SD 的 HR ================\n")
print(out, row.names = FALSE, digits = 3)
write.table(out, paste0(OUT, "/output/survival_metric_battery.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n================ 判定 ================\n")
for (ep in c("OS", "PFI")) {
  s <- out[out$endpoint == ep, ]
  sig <- s[s$p_adj_within_endpoint < 0.05, ]
  cat(sprintf("%s: %d 个指标中，BH<0.05 者 %d 个\n", ep, nrow(s), nrow(sig)))
  if (nrow(sig)) print(sig[, c("metric", "HR_per_SD", "p", "p_adj_within_endpoint")],
                       row.names = FALSE, digits = 3)
}

## 多变量 Cox（修复 NA / 单水平因子问题）
cat("\n================ 多变量 Cox（校正年龄/性别/亚型）================\n")
d <- tab
d$Subtype <- droplevels(factor(d$Subtype))
d$sex <- droplevels(factor(d$sex))
for (ep in c("OS", "PFI")) {
  dd <- d[!is.na(d[[ep]]) & !is.na(d[[paste0(ep, ".time")]]) &
          !is.na(d$burden_sd) & !is.na(d$age) & !is.na(d$sex) &
          !is.na(d$Subtype), ]
  dd$Subtype <- droplevels(dd$Subtype); dd$sex <- droplevels(dd$sex)
  cat(sprintf("  %s: 完整病例 n=%d, 亚型水平=%s\n", ep, nrow(dd),
              paste(levels(dd$Subtype), collapse = "/")))
  m <- coxph(as.formula(paste0("Surv(", ep, ".time, ", ep,
                               ") ~ burden_sd + age + sex + Subtype")), data = dd)
  s <- summary(m)
  print(round(cbind(s$coefficients[, c("exp(coef)", "Pr(>|z|)")],
                    s$conf.int[, c("lower .95", "upper .95")]), 4))
}
