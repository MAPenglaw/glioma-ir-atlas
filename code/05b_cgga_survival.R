#!/usr/bin/env Rscript
# Step 5b — CGGA 外部验证：RI 轴表达签名的临床意义（n=970）

suppressMessages({ library(survival); library(data.table) })

OUT <- "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
cg <- read.delim(paste0(OUT, "/output/cgga_riaxis_scores.tsv"),
                 row.names = 1, check.names = FALSE)
cg$RIaxis_score <- as.numeric(cg$RIaxis_score)
cg$z <- (cg$RIaxis_score - mean(cg$RIaxis_score)) / sd(cg$RIaxis_score)
cg$OSy <- cg$OS_months / 12
cg$grade <- factor(cg$Grade, levels = c("II", "III", "IV"))
cg$idh <- factor(cg$IDH)
cg$sex <- factor(cg$Gender)
cg$grp <- factor(ifelse(cg$z > median(cg$z), "RIaxis-high", "RIaxis-low"),
                 levels = c("RIaxis-low", "RIaxis-high"))
cg$codel <- factor(cg$codel_1p19q)

cat("CGGA 队列: n =", nrow(cg), ", 事件 =", sum(cg$event), "\n")
cat("随访中位（年）:", round(median(cg$OSy), 2), "\n")

## ---------- 1. 签名与临床病理特征的关系 ----------
cat("\n=========== RI 轴签名 vs 临床病理特征 ===========\n")
for (v in c("idh", "grade", "sex", "codel")) {
  d <- cg[!is.na(cg[[v]]), ]
  if (nlevels(droplevels(d[[v]])) < 2) next
  k <- kruskal.test(d$RIaxis_score ~ droplevels(d[[v]]))
  cat(sprintf("%-8s Kruskal-Wallis p = %.3e\n", v, k$p.value))
  print(tapply(d$RIaxis_score, droplevels(d[[v]]), function(x)
    round(c(n = length(x), mean = mean(x)), 3)))
}

## ---------- 2. 总体生存 ----------
cat("\n=========== 单变量 Cox（整队列）===========\n")
uni <- function(d, label) {
  if (nrow(d) < 40) { cat(sprintf("%-28s n 不足 (%d)\n", label, nrow(d))); return(NULL) }
  m <- coxph(Surv(OSy, event) ~ z, data = d)
  s <- summary(m)
  cat(sprintf("%-28s n=%4d ev=%3d HR/SD=%.3f [%.3f-%.3f] p=%.3e\n", label,
              m$n, m$nevent, s$coefficients[1, "exp(coef)"],
              s$conf.int[1, "lower .95"], s$conf.int[1, "upper .95"],
              s$coefficients[1, "Pr(>|z|)"]))
  data.frame(subset = label, n = m$n, events = m$nevent,
             HR_per_SD = s$coefficients[1, "exp(coef)"],
             lower = s$conf.int[1, "lower .95"],
             upper = s$conf.int[1, "upper .95"],
             p = s$coefficients[1, "Pr(>|z|)"])
}
res <- uni(cg, "全队列")
res <- rbind(res, uni(cg[cg$grade == "IV", ], "WHO IV (GBM)")
                 , uni(cg[cg$grade == "III", ], "WHO III")
                 , uni(cg[cg$grade == "II", ], "WHO II")
                 , uni(cg[cg$idh == "Wildtype", ], "IDH-wildtype")
                 , uni(cg[cg$idh == "Mutant", ], "IDH-mutant")
                 , uni(cg[cg$idh == "Wildtype" & cg$grade == "IV", ], "IDH-wt WHO IV")
                 , uni(cg[cg$codel == "Codel", ], "1p19q codel")
                 , uni(cg[cg$codel == "Non-codel", ], "1p19q non-codel"))

## ---------- 3. Log-rank（中位分割）----------
cat("\n=========== Log-rank（中位分割）===========\n")
lr <- function(d, label) {
  if (nrow(d) < 40) return(NULL)
  s <- survdiff(Surv(OSy, event) ~ grp, data = d)
  p <- 1 - pchisq(s$chisq, 1)
  cat(sprintf("%-28s chisq=%.3f p=%.4f\n", label, s$chisq, p))
  data.frame(subset = label, n = nrow(d), chisq = s$chisq, p = p)
}
lrres <- rbind(lr(cg, "全队列"),
               lr(cg[cg$grade == "IV", ], "WHO IV (GBM)"),
               lr(cg[cg$idh == "Wildtype", ], "IDH-wildtype"),
               lr(cg[cg$idh == "Mutant", ], "IDH-mutant"))

## ---------- 4. 多变量 Cox（校正分级/IDH/年龄/性别/1p19q）----------
cat("\n=========== 多变量 Cox ===========\n")
run_mv <- function(d, formula, label) {
  d <- droplevels(d)
  m <- coxph(formula, data = d)
  s <- summary(m)
  cat("\n--", label, " n =", m$n, " events =", m$nevent, "\n")
  print(round(cbind(s$coefficients[, c("exp(coef)", "Pr(>|z|)")],
                    s$conf.int[, c("lower .95", "upper .95")]), 4))
  data.frame(subset = label, term = rownames(s$coefficients),
             HR = s$coefficients[, "exp(coef)"],
             lower = s$conf.int[, "lower .95"],
             upper = s$conf.int[, "upper .95"],
             p = s$coefficients[, "Pr(>|z|)"], n = m$n, events = m$nevent)
}
d1 <- cg[!is.na(cg$grade) & !is.na(cg$idh) & !is.na(cg$Age) & !is.na(cg$sex), ]
mv1 <- run_mv(d1, Surv(OSy, event) ~ z + grade + idh + Age + sex,
              "全队列 校正分级/IDH/年龄/性别")
d2 <- d1[!is.na(d1$codel), ]
mv2 <- run_mv(d2, Surv(OSy, event) ~ z + grade + idh + Age + sex + codel,
              "全队列 +1p19q")
d3 <- d1[d1$idh == "Wildtype" & d1$grade == "IV", ]
mv3 <- run_mv(d3, Surv(OSy, event) ~ z + Age + sex, "IDH-wt WHO IV")

## ---------- 5. 中位生存 ----------
cat("\n=========== 中位生存（月）===========\n")
for (lab in c("全队列", "WHO IV (GBM)", "IDH-wildtype")) {
  d <- switch(lab, "全队列" = cg, "WHO IV (GBM)" = cg[cg$grade == "IV", ],
              "IDH-wildtype" = cg[cg$idh == "Wildtype", ])
  s <- survfit(Surv(OSy, event) ~ grp, data = d)
  cat("--", lab, "\n"); print(s)
}

write.table(res, paste0(OUT, "/output/cgga_cox_univariate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(mv1, mv2, mv3), paste0(OUT, "/output/cgga_cox_multivariate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lrres, paste0(OUT, "/output/cgga_logrank.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cg, paste0(OUT, "/output/cgga_analysis.rds"))
cat("\n完成 -> cgga_cox_univariate.tsv / cgga_cox_multivariate.tsv / cgga_logrank.tsv\n")
