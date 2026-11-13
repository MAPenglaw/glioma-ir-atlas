#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G2 — 分基因层面：剪接因子表达 vs 全局 RI 负荷

把 ssGSEA 的「基因集层面」结论下沉到「单个基因」层面，
重点回答：SNRPC 及 U1 snRNP 各组分与 RI 负荷的关系是什么？

设计：
  1. 全剪接因子池（214 个，与方案 A 同源）逐个 Spearman + BH
  2. U1 snRNP panel（12 个，与 03b ssGSEA 完全同源）与核心剪接体 panel（12 个）定向检验
  3. 组内一致性检验：U1 panel 内负相关比例是否 > 50%（二项 + 置换零分布）
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
BASE = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = f"{BASE}/output"
POOL = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/sf_signature/output/sf_gene_pool.json"
N_PERM = 20000
SEED = 20260917

U1_SNRNP = ["SNRPC", "SNRPA", "SNRPB", "SNRPD1", "SNRPD2", "SNRPD3",
            "SNRPE", "SNRPF", "SNRPG", "SNRNP70", "SNRNP200", "PRPF8"]
SPLICE_CORE = ["SF3B1", "SF3B2", "SF3A1", "U2AF1", "U2AF2", "SRSF1",
               "SRSF2", "SRSF3", "HNRNPA1", "HNRNPC", "PTBP1", "RBM39"]

# GBM_HiSeqV2（hg19/Xena）用的是旧 HGNC 符号，SR 蛋白等基因会静默缺失
ALIAS = {
    "SRSF1": "SFRS1", "SRSF2": "SFRS2", "SRSF3": "SFRS3", "SRSF4": "SFRS4",
    "SRSF5": "SFRS5", "SRSF6": "SFRS6", "SRSF7": "SFRS7", "SRSF8": "SFRS8",
    "SRSF9": "SFRS9", "SRSF10": "SFRS13A", "SRSF11": "SFRS11", "SRSF12": "SFRS12",
    "SUGP2": "SFRS14", "SCAF4": "SFRS15", "SCAF11": "SFRS16", "AKAP17A": "SFRS17A",
    "CPSF5": "NUDT21", "RBFOX1": "A2BP1", "RBFOX2": "RBM9",
}


def barcode(c):
    return str(c).strip().strip('"')[:15]


def bh(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def main():
    print("=" * 74)
    print("G2  剪接因子表达 vs 全局 RI 负荷（分基因层面）")
    print("=" * 74)

    ann = pd.read_csv(f"{OUT}/sample_RI_group_annotated.tsv", sep="\t", index_col=0)
    ann.index = [str(i)[:15] for i in ann.index]
    burden = ann["RI_burden"]
    print(f"[1] RI 负荷样本: {len(burden)} 例")

    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0, compression="gzip")
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[~expr.index.duplicated()]
    common = [c for c in burden.index if c in set(expr.columns)]
    E = expr[common]
    b = burden.loc[common]
    print(f"[2] 表达与负荷交集样本: {len(common)} 例")

    pool = json.load(open(POOL))["pool"]
    panel = sorted(set(pool) | set(U1_SNRNP) | set(SPLICE_CORE))
    gi = set(E.index)
    pairs = []          # (标准名, 矩阵中实际行名)
    for g in panel:
        if g in gi:
            pairs.append((g, g))
        elif ALIAS.get(g) in gi:
            pairs.append((g, ALIAS[g]))
    miss = [g for g in panel if g not in gi and ALIAS.get(g) not in gi]
    aliased = [(g, a) for g, a in pairs if g != a]
    print(f"[3] 候选剪接因子 {len(panel)} 个｜可测 {len(pairs)} 个｜缺失 {len(miss)} 个")
    if aliased:
        print(f"    经旧名映射命中 {len(aliased)} 个：" +
              ", ".join(f"{g}→{a}" for g, a in aliased))
    if miss:
        print("    仍缺失:", ", ".join(miss))

    rows = []
    for g, row in pairs:
        v = E.loc[row].astype(float)
        m = v.notna() & b.notna()
        r = stats.spearmanr(v[m], b[m])
        # 组间
        hi = b[m] > b[m].median()
        mw = stats.mannwhitneyu(v[m][hi], v[m][~hi], alternative="two-sided")
        rows.append({"gene": g, "matrix_row": row, "rho": r.statistic, "p": r.pvalue,
                     "delta_high_low": v[m][hi].mean() - v[m][~hi].mean(),
                     "p_mw": mw.pvalue})
    df = pd.DataFrame(rows)
    df["q"] = bh(df["p"].values)
    df = df.sort_values("rho")
    df.to_csv(f"{OUT}/G2_sf_vs_burden.tsv", sep="\t", index=False)

    n_sig = int((df["q"] < 0.05).sum())
    print(f"\n[4] 全池结果：{len(df)} 个基因｜BH q<0.05: {n_sig} 个")
    print("    --- 与 RI 负荷负相关最强的前 12 个 ---")
    print(df.head(12)[["gene", "rho", "p", "q"]].to_string(index=False,
          float_format=lambda x: f"{x:.4g}"))
    print("    --- 正相关最强的前 8 个 ---")
    print(df.tail(8)[["gene", "rho", "p", "q"]].iloc[::-1].to_string(index=False,
          float_format=lambda x: f"{x:.4g}"))

    # ---------- U1 snRNP panel 定向检验 ----------
    print("\n[5] U1 snRNP panel（与 03b ssGSEA 同源）逐个检验")
    sub = df[df["gene"].isin(U1_SNRNP)].copy()
    print(sub[["gene", "rho", "p", "q", "delta_high_low"]].to_string(index=False,
          float_format=lambda x: f"{x:.4g}"))
    n_neg = int((sub["rho"] < 0).sum())
    bt = stats.binomtest(n_neg, len(sub), 0.5, alternative="two-sided")
    print(f"    负相关 {n_neg}/{len(sub)}｜二项检验 p = {bt.pvalue:.4g}"
          f"｜中位 ρ = {sub['rho'].median():+.4f}")

    # 置换零分布：随机抽同样大小的基因集，看负相关比例
    rng = np.random.default_rng(SEED)
    all_rho = df["rho"].values
    null_frac = np.array([np.mean(rng.choice(all_rho, len(sub), replace=False) < 0)
                          for _ in range(N_PERM)])
    p_perm = (np.sum(np.abs(null_frac - 0.5) >= abs(n_neg / len(sub) - 0.5)) + 1) / (N_PERM + 1)
    print(f"    置换零分布（从全池随机抽 {len(sub)} 个，{N_PERM} 次）："
          f"负相关比例 {null_frac.mean():.3f} ± {null_frac.std():.3f} → p = {p_perm:.4g}")

    # ---------- 核心剪接体 panel ----------
    print("\n[6] 核心剪接体 panel（预期代偿性上调 → 正相关）")
    sub2 = df[df["gene"].isin(SPLICE_CORE)].copy()
    print(sub2[["gene", "rho", "p", "q", "delta_high_low"]].to_string(index=False,
          float_format=lambda x: f"{x:.4g}"))
    n_pos = int((sub2["rho"] > 0).sum())
    bt2 = stats.binomtest(n_pos, len(sub2), 0.5, alternative="two-sided")
    print(f"    正相关 {n_pos}/{len(sub2)}｜二项检验 p = {bt2.pvalue:.4g}"
          f"｜中位 ρ = {sub2['rho'].median():+.4f}")

    # ---------- SNRPC 单独 ----------
    print("\n[7] SNRPC 单独检验")
    s = df[df["gene"] == "SNRPC"]
    if len(s):
        s = s.iloc[0]
        print(f"    ρ = {s['rho']:+.4f}｜p = {s['p']:.4g}｜BH q = {s['q']:.4g}"
              f"｜Δ(RI-high − RI-low) = {s['delta_high_low']:+.4f}")
        print(f"    全池排名：按 ρ 升序第 {int((df['rho'] < s['rho']).sum()) + 1}/{len(df)} 位（越靠前越负）")

    print(f"\n[8] 已输出 {OUT}/G2_sf_vs_burden.tsv")


if __name__ == "__main__":
    main()
