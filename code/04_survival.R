#!/usr/bin/env Rscript
# Step 4 — RI 负荷 / RI 分型的生存分析（OS / DSS / PFI）
# 这是决定整个设计成立与否的一步：若 RI 轴无预后价值，则「预后」框架不成立。

suppressMessages({ library(survival); library(data.table) })

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"

tab <- read.delim(paste0(OUT, "/output/sample_RI_group_clinical.tsv"),
                  row.names = 1, check.names = FALSE)
gs  <- read.delim(paste0(OUT, "/output/gsva_scores.tsv"),
                  row.names = 1, check.names = FALSE)
tab <- cbind(tab, gs[rownames(tab), c("U1_snRNP", "Spliceosome_core",
                                      "Proliferation", "Glioma_StemCell",
                                      "IFN_gamma_signature", "Treg", "CD8_T_cell")])

## 亚型编码反推结果（marker 验证）: 0=Classical 1=Neural 2=Proneural 3=Mesenchymal
SUBT <- c("0" = "Classical", "1" = "Neural", "2" = "Proneural", "3" = "Mesenchymal")
tab$Subtype <- SUBT[as.character(tab$GeneExp_Subtype)]
tab$age <- suppressWarnings(as.numeric(tab$age_at_initial_pathologic_diagnosis))
## 性别编码反推结果（XIST/KDM5D/RPS4Y1 验证）: 0=Female, 2=Male
tab$sex <- factor(ifelse(tab$gender == 2, "Male",
                  ifelse(tab$gender == 0, "Female", NA)),
                  levels = c("Female", "Male"))
tab$RI_group <- factor(tab$RI_group, levels = c("RI-low", "RI-high"))
tab$RI_burden <- as.numeric(tab$RI_burden)
tab$burden_sd <- (tab$RI_burden - mean(tab$RI_burden, na.rm = TRUE)) /
                 sd(tab$RI_burden, na.rm = TRUE)
tab$burden_tertile <- cut(tab$RI_burden,
                          breaks = quantile(tab$RI_burden, c(0, 1/3, 2/3, 1)),
                          labels = c("T1(low)", "T2", "T3(high)"),
                          include.lowest = TRUE)
cat("样本数:", nrow(tab), " 有 OS:", sum(!is.na(tab$OS)), "\n")
print(table(tab$RI_group, tab$Subtype, useNA = "ifany"))

cat("\n============ 单变量 Cox（连续 RI 负荷，每 +1 SD）============\n")
cox_tab <- data.frame()
for (ep in c("OS", "DSS", "PFI")) {
  for (v in c("burden_sd", "RI_group", "age", "sex", "Subtype",
              "U1_snRNP", "Spliceosome_core", "Proliferation",
              "Glioma_StemCell", "IFN_gamma_signature")) {
    d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]) &
             !is.na(tab[[v]]), ]
    if (v == "age") d <- d[!is.na(d$age), ]
    if (nrow(d) < 40) next
    f <- as.formula(paste0("Surv(", ep, ".time, ", ep, ") ~ ", v))
    m <- tryCatch(coxph(f, data = d), error = function(e) NULL)
    if (is.null(m)) next
    s <- summary(m)
    co <- s$coefficients; ci <- s$conf.int
    # 多水平变量只取第一行做汇总
    cox_tab <- rbind(cox_tab, data.frame(
      endpoint = ep, variable = v,
      n = m$n, events = m$nevent,
      HR = co[1, "exp(coef)"], lower = ci[1, "lower .95"],
      upper = ci[1, "upper .95"], p = co[1, "Pr(>|z|)"]))
  }
}
cox_tab$p_adj <- p.adjust(cox_tab$p, method = "BH")
print(cox_tab[cox_tab$variable %in% c("burden_sd", "RI_group"), ],
      row.names = FALSE, digits = 4)

cat("\n============ 多变量 Cox（校正年龄/性别/分子亚型）============\n")
mv_tab <- data.frame()
for (ep in c("OS", "DSS", "PFI")) {
  d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]) &
           !is.na(tab$burden_sd) & !is.na(tab$age) & !is.na(tab$sex) &
           !is.na(tab$Subtype), ]
  d$Subtype <- droplevels(factor(d$Subtype))
  d$sex <- droplevels(factor(d$sex))
  if (nrow(d) < 40 || nlevels(d$Subtype) < 2) {
    cat(sprintf("  %s: 完整病例不足或亚型单一（n=%d），跳过\n", ep, nrow(d)))
    next
  }
  m <- coxph(as.formula(paste0("Surv(", ep, ".time, ", ep,
                               ") ~ burden_sd + age + sex + Subtype")), data = d)
  s <- summary(m)
  for (i in seq_len(nrow(s$coefficients))) {
    mv_tab <- rbind(mv_tab, data.frame(
      endpoint = ep, term = rownames(s$coefficients)[i],
      n = m$n, events = m$nevent,
      HR = s$coefficients[i, "exp(coef)"],
      lower = s$conf.int[i, "lower .95"],
      upper = s$conf.int[i, "upper .95"],
      p = s$coefficients[i, "Pr(>|z|)"]))
  }
}
print(mv_tab[mv_tab$term == "burden_sd", ], row.names = FALSE, digits = 4)

cat("\n============ Log-rank（RI 分组 / 三分位）============\n")
lr <- data.frame()
for (ep in c("OS", "DSS", "PFI")) {
  for (v in c("RI_group", "burden_tertile")) {
    d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]) &
             !is.na(tab[[v]]), ]
    sd <- survdiff(as.formula(paste0("Surv(", ep, ".time, ", ep, ") ~ ", v)), data = d)
    p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    cat(sprintf("  %s ~ %-16s n=%d chisq=%.3f p=%.4f\n", ep, v, nrow(d), sd$chisq, p))
    lr <- rbind(lr, data.frame(endpoint = ep, variable = v, n = nrow(d),
                               chisq = sd$chisq, p = p))
  }
}

cat("\n============ 中位生存（RI 分组）============\n")
for (ep in c("OS", "DSS", "PFI")) {
  d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]), ]
  s <- survfit(as.formula(paste0("Surv(", ep, ".time, ", ep,
                                 ") ~ RI_group")), data = d)
  print(s)
}

cat("\n============ 连续 RI 负荷分层（三分位）KM 中位值 ============\n")
for (ep in c("OS", "PFI")) {
  d <- tab[!is.na(tab[[ep]]) & !is.na(tab[[paste0(ep, ".time")]]), ]
  s <- survfit(as.formula(paste0("Surv(", ep, ".time, ", ep,
                                 ") ~ burden_tertile")), data = d)
  print(s)
}

write.table(cox_tab, paste0(OUT, "/output/survival_cox_univariate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mv_tab, paste0(OUT, "/output/survival_cox_multivariate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr, paste0(OUT, "/output/survival_logrank.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(tab, paste0(OUT, "/output/survival_input.rds"))
cat("\n完成 -> survival_cox_univariate.tsv / survival_cox_multivariate.tsv / survival_logrank.tsv\n")
