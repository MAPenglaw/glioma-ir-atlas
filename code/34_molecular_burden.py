#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""34 — B-2：IR 轴与分子分型的关联。

把 B 线已有的 IR 轴（output_km/{GBM,LGG}/pc1_burden.tsv：PC1 与 RI burden）
与 33 号获得的分子分型（IDH1/2 突变、CIC/FUBP1 共缺失代理）对齐，回答：
「这条既不是预后、也不是增殖替身的轴，到底对应什么生物学？」

严格遵守项目红线：
  - 绝不跨 GBM/LGG 混合比较（R3/R4），一律在队列内分层
  - 不谈预后，只谈横断面的组间差异
  - 报告效应量（Cliff's delta 与 bootstrap 95% CI），不以 p 值论有无

输出：output_pancan/mol_burden_stats.tsv / mol_burden_plot.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
KM = os.path.join(ROOT, "output_km")
OUT = os.path.join(ROOT, "output_pancan")
CLIN = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/data/TCGA_GBMLGG_clinical.tsv"

rng = np.random.default_rng(20260928)


def cliffs_delta(a, b, n_boot=2000):
    """Cliff's delta = 2*AUC - 1，a 组相对 b 组；返回 delta 与 bootstrap 95% CI。"""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan, (np.nan, np.nan)
    U = mannwhitneyu(a, b, alternative="two-sided").statistic
    auc = U / (n1 * n2)
    delta = 2 * auc - 1
    boots = []
    for _ in range(n_boot):
        aa = rng.choice(a, n1, replace=True)
        bb = rng.choice(b, n2, replace=True)
        Ub = mannwhitneyu(aa, bb, alternative="two-sided").statistic
        boots.append(2 * Ub / (n1 * n2) - 1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(delta), (float(lo), float(hi))


def main():
    mol = pd.read_csv(os.path.join(OUT, "glio_molecular.tsv"), sep="\t")
    frames = []
    for c in ("GBM", "LGG"):
        pb = pd.read_csv(os.path.join(KM, c, "pc1_burden.tsv"), sep="\t")
        pb["cohort"] = c
        frames.append(pb)
    ir = pd.concat(frames)

    d = ir.merge(mol, left_on=["sample", "cohort"],
                 right_on=["sample_id", None], how="left") \
        if False else ir.merge(
            mol.rename(columns={"sample_id": "sample"}), on="sample", how="left")

    # 补 WHO 分级（本地临床表）
    clin = pd.read_csv(CLIN, sep="\t", dtype=str, low_memory=False)
    clin["patient"] = clin["sampleID"].astype(str).str[:12]
    gcol = "neoplasm_histologic_grade"
    if gcol in clin.columns:
        gr = (clin[["patient", gcol]].dropna()
              .drop_duplicates("patient"))
        gr.columns = ["patient", "grade"]
        d["patient"] = d["sample"].str[:12]
        d = d.merge(gr, on="patient", how="left")
    else:
        d["grade"] = np.nan

    d = d.dropna(subset=["IDH_status"])
    print(f"[merge] 有 IR 轴且分子分型齐全的样本 {len(d)}", flush=True)
    print(d.groupby(["cohort", "IDH_status"]).size().to_string(), flush=True)

    rows = []
    for c in ("GBM", "LGG"):
        sub = d[d.cohort == c]
        if sub.IDH_status.nunique() < 2:
            continue
        for metric in ("RI_burden", "PC1"):
            a = sub.loc[sub.IDH_status == "IDHmut", metric].dropna()
            b = sub.loc[sub.IDH_status == "IDHwt", metric].dropna()
            if len(a) < 3 or len(b) < 3:
                continue
            U = mannwhitneyu(a, b, alternative="two-sided")
            delta, (lo, hi) = cliffs_delta(a, b)
            rows.append({
                "cohort": c, "metric": metric, "comparison": "IDHmut vs IDHwt",
                "n1": len(a), "n2": len(b),
                "median1": round(float(a.median()), 4),
                "median2": round(float(b.median()), 4),
                "cliffs_delta": round(delta, 3),
                "ci_lo": round(lo, 3), "ci_hi": round(hi, 3),
                "p": float(U.pvalue),
            })
            print(f"[{c}] {metric}: IDHmut n={len(a)} 中位 {a.median():.4f} | "
                  f"IDHwt n={len(b)} 中位 {b.median():.4f} | "
                  f"delta={delta:.3f} [{lo:.3f},{hi:.3f}] p={U.pvalue:.3g}",
                  flush=True)

        # LGG 内三分类
        if c == "LGG":
            for metric in ("RI_burden", "PC1"):
                groups = [g[metric].dropna().values
                          for _, g in sub.groupby("mol3")]
                labels = [k for k, _ in sub.groupby("mol3")]
                if len(groups) < 2:
                    continue
                H = kruskal(*groups)
                print(f"[LGG] {metric} Kruskal-Wallis 三分类 "
                      f"H={H.statistic:.2f} p={H.pvalue:.3g}", flush=True)
                print("     各组中位数: " + ", ".join(
                    f"{l}={np.median(g):.4f}(n={len(g)})"
                    for l, g in zip(labels, groups)), flush=True)
                rows.append({
                    "cohort": c, "metric": metric,
                    "comparison": "Kruskal 3-class (codel/noncodel/IDHwt)",
                    "n1": len(groups), "n2": np.nan,
                    "median1": np.nan, "median2": np.nan,
                    "cliffs_delta": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                    "p": float(H.pvalue),
                })
                # 两两
                for i in range(len(labels)):
                    for j in range(i + 1, len(labels)):
                        dl, (l2, h2) = cliffs_delta(groups[i], groups[j])
                        U = mannwhitneyu(groups[i], groups[j],
                                         alternative="two-sided")
                        rows.append({
                            "cohort": c, "metric": metric,
                            "comparison": f"{labels[i]} vs {labels[j]}",
                            "n1": len(groups[i]), "n2": len(groups[j]),
                            "median1": round(float(np.median(groups[i])), 4),
                            "median2": round(float(np.median(groups[j])), 4),
                            "cliffs_delta": round(dl, 3),
                            "ci_lo": round(l2, 3), "ci_hi": round(h2, 3),
                            "p": float(U.pvalue),
                        })
                        print(f"      {labels[i]} vs {labels[j]}: "
                              f"delta={dl:.3f} [{l2:.3f},{h2:.3f}] "
                              f"p={U.pvalue:.3g}", flush=True)

    # LGG 内按 WHO 分层的 IDH 效应（回答：是 IDH 的替身还是分级的替身）
    if "grade" in d.columns:
        print("\n=== LGG 内按 WHO 分层的 IDHmut vs IDHwt ===", flush=True)
        sub = d[(d.cohort == "LGG") & d.grade.notna()]
        for gv, gg in sub.groupby("grade"):
            a = gg.loc[gg.IDH_status == "IDHmut", "RI_burden"].dropna()
            b = gg.loc[gg.IDH_status == "IDHwt", "RI_burden"].dropna()
            if len(a) < 3 or len(b) < 3:
                print(f"  grade={gv}: IDHmut n={len(a)}, IDHwt n={len(b)} "
                      f"样本不足，跳过", flush=True)
                continue
            U = mannwhitneyu(a, b, alternative="two-sided")
            dl, (l2, h2) = cliffs_delta(a, b)
            rows.append({
                "cohort": "LGG", "metric": "RI_burden",
                "comparison": f"IDHmut vs IDHwt (within grade {gv})",
                "n1": len(a), "n2": len(b),
                "median1": round(float(a.median()), 4),
                "median2": round(float(b.median()), 4),
                "cliffs_delta": round(dl, 3),
                "ci_lo": round(l2, 3), "ci_hi": round(h2, 3),
                "p": float(U.pvalue),
            })
            print(f"  grade={gv}: IDHmut n={len(a)} 中位 {a.median():.4f} | "
                  f"IDHwt n={len(b)} 中位 {b.median():.4f} | "
                  f"delta={dl:.3f} [{l2:.3f},{h2:.3f}] p={U.pvalue:.3g}",
                  flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "mol_burden_stats.tsv"), sep="\t", index=False)
    keep = [c for c in ["sample", "cohort", "RI_burden", "PC1", "IDH_status",
                        "mol3", "codel_proxy", "grade", "IDH1", "IDH2",
                        "TP53", "ATRX", "EGFR", "CIC", "FUBP1"]
            if c in d.columns]
    d[keep].to_csv(os.path.join(OUT, "mol_burden_plot.tsv"), sep="\t",
                   index=False)
    print(f"\n[sav] mol_burden_stats.tsv ({len(df)} 行) / "
          f"mol_burden_plot.tsv ({len(d)} 行)", flush=True)


if __name__ == "__main__":
    main()
