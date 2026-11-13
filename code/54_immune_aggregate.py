#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""54 — 免疫综合评分（Immune_aggregate）作为主结果（评审 P1-3）。

问题：17 个免疫/基质面板高度相关，"17/17 显著"不是 17 次独立检验，
     不能作为效应量证据。做法：把 15 个免疫面板（排除基质 CAF_Stromal
     与内皮 Endothelial 两个非免疫谱系）的 z 分数取均值，得到单一
     Immune_aggregate，作为**主结果**；17 个面板退为一致性佐证。

输出：output_pancan/immune_aggregate.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
OUTP = os.path.join(ROOT, "ri_atlas", "output_pancan")

IMMUNE = ["T_cell_general", "CD8_T", "Treg", "NK", "B_cell", "Myeloid_general",
          "Microglia", "M2_like", "Neutrophil", "DC", "Cytotoxic",
          "IFN_gamma_Ayers", "MHC_I", "MHC_II", "Checkpoint"]
NON_IMMUNE = ["CAF_Stromal", "Endothelial"]

MODELS = [
    ("M0_raw", []),
    ("M1_TotalExpr", ["TotalExpr"]),
    ("M2_plus_Neuronal", ["TotalExpr", "Neuronal_normal_brain"]),
    ("M3_plus_Spliceosome", ["TotalExpr", "Neuronal_normal_brain", "Sm_core",
                             "Proliferation"]),
    ("M4_plus_Purity", ["TotalExpr", "Neuronal_normal_brain", "Sm_core",
                        "Proliferation", "Tumor_content"]),
    ("M5_plus_IDH", ["TotalExpr", "Neuronal_normal_brain", "Sm_core",
                     "Proliferation", "Tumor_content", "IDHmut"]),
    ("M7_hi_lo_layers", ["TotalExpr_hi", "TotalExpr_lo", "Neuronal_normal_brain",
                         "Sm_core", "Proliferation", "Tumor_content", "IDHmut"]),
]


def ranknorm(v):
    v = np.asarray(v, dtype=float)
    return st.norm.ppf(st.rankdata(v) / (len(v) + 1.0))


def partial_corr(x, y, covs):
    """秩偏相关：rank→正态分位后 OLS 残差 Pearson。k=0 时即 Spearman。"""
    xr, yr = ranknorm(x), ranknorm(y)
    if len(covs) == 0:
        r, p = st.pearsonr(xr, yr)
        n_eff = len(x)
    else:
        C = np.column_stack([ranknorm(c) for c in covs] + [np.ones(len(x))])
        bx, *_ = np.linalg.lstsq(C, xr, rcond=None)
        by, *_ = np.linalg.lstsq(C, yr, rcond=None)
        rx, ry = xr - C @ bx, yr - C @ by
        r, p = st.pearsonr(rx, ry)
        n_eff = len(x) - C.shape[1] + 1
    # Fisher z 95% CI
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1 / np.sqrt(max(n_eff - 3, 1))
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    return float(r), float(lo), float(hi), float(p)


def main():
    sc = pd.read_csv(os.path.join(OUTP, "tme_sample_scores.tsv"), sep="\t")
    cf = pd.read_csv(os.path.join(OUTP, "tme_confounders.tsv"), sep="\t")
    cf["key"] = cf["sample"].astype(str).str[:12]

    rows = []
    for cohort in ("GBM", "LGG"):
        g = sc[(sc["cohort"] == cohort) & (sc["panel"].isin(IMMUNE))]
        imm = [p for p in IMMUNE if p in set(g["panel"])]
        absent = [p for p in IMMUNE if p not in set(g["panel"])]
        if absent:
            print(f"   [note] 面板缺失，改用 {len(imm)} 个：缺 {absent}", flush=True)
        wide = g.pivot_table(index="sample", columns="panel", values="score",
                             aggfunc="first")
        # 面板内部一致性（Cronbach alpha），说明"17 次检验"确实不独立
        X = wide[imm].dropna()
        k = X.shape[1]
        alpha = k / (k - 1) * (1 - X.var(axis=0, ddof=1).sum() /
                               X.sum(axis=1).var(ddof=1))
        agg = wide[imm].mean(axis=1).rename("Immune_aggregate")

        bur = g.groupby("sample")["RI_burden"].first()
        df = pd.concat([agg, bur], axis=1).dropna()
        df["key"] = df.index.astype(str).str[:12]
        # 混杂表中的 RI_burden 与面板表重复，排除以免列名冲突
        cf2 = cf.drop(columns=[c for c in ("RI_burden", "cohort", "sample")
                               if c in cf.columns])
        df = df.merge(cf2.drop_duplicates("key").set_index("key"),
                      left_on="key", right_index=True, how="left")
        df = df.dropna(subset=["RI_burden", "Immune_aggregate", "TotalExpr"])
        print(f"[{cohort}] n={len(df)}；{len(imm)} 面板 Cronbach α = {alpha:.3f}",
              flush=True)

        for mname, covs in MODELS:
            use = [c for c in covs if c in df.columns]
            miss = [c for c in covs if c not in df.columns]
            if miss:
                print(f"   [skip] {mname} 缺 {miss}", flush=True)
                continue
            ok = df[use + ["RI_burden", "Immune_aggregate"]].notna().all(axis=1)
            sub = df[ok]
            r, lo, hi, p = partial_corr(sub["RI_burden"].values,
                                        sub["Immune_aggregate"].values,
                                        [sub[c].values for c in use])
            rows.append({"cohort": cohort, "panel": "Immune_aggregate",
                         "model": mname, "n": len(sub), "rho": r,
                         "ci_lo": lo, "ci_hi": hi, "p": p,
                         "k_cov": len(use),
                         "cronbach_alpha": round(float(alpha), 4)})
            print(f"   {mname:<20} n={len(sub):3d} rho={r:+.3f} "
                  f"({lo:+.3f},{hi:+.3f}) p={p:.2e}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUTP, "immune_aggregate.tsv"), sep="\t", index=False)
    print("DONE -> output_pancan/immune_aggregate.tsv")


if __name__ == "__main__":
    main()
