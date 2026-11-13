#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 关键对照 —— 这个负相关是生物学（NMD）还是 reads 分配的数学耦合？

问题：IR 事件的 PSI 基于「内含子 reads / (内含子 + 外显子 reads)」，
      而基因表达定量主要来自外显子 reads。内含子保留升高会把 reads 从外显子
      「挪」到内含子，从而机械性地压低基因定量 —— 这与 NMD 的预期方向完全一致，
      两者无法从 IR 事件本身区分。

判别设计（阴性对照）：
      对**不涉及内含子**的剪接事件类型做完全相同的分析。
        · ME（互斥外显子）—— 纯外显子之间的选择，定量上不应有 reads 转移
        · ES（外显子跳跃）—— 同上，只是外显子的进/出
        · A3 / A5（可变 3'/5' 位点）—— 边界微调
      若这些事件同样显示 ~65% 负相关，则 IR 的负相关主要是技术耦合；
      若只有 IR 显著偏离 50%，则支持内含子滞留特异的 NMD 解释。

附加对照：IR 事件按「内含子长度代理 = 检出 PSI 均值」分层，
         高 PSI（滞留多、reads 转移多）事件若负相关更强 → 支持技术耦合。
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
BASE = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = f"{BASE}/output"

MIN_DETECT = 0.90
N_PERM = 2000
SEED = 20260917


def barcode(c):
    return str(c).strip().strip('"')[:15]


def rank_rows(M):
    return M.apply(lambda r: r.rank(), axis=1)


def row_corr(A, B):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    a = np.nan_to_num(A - np.nanmean(A, axis=1, keepdims=True))
    b = np.nan_to_num(B - np.nanmean(B, axis=1, keepdims=True))
    num = (a * b).sum(axis=1)
    den = np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))
    return num / np.where(den == 0, np.nan, den)


def analyze(Pk, Ek, label, n_perm=N_PERM, seed=SEED):
    P_fill = Pk.T.fillna(Pk.median(axis=1)).T
    Rp = rank_rows(P_fill).values
    Re = rank_rows(Ek).values
    rho = row_corr(Rp, Re)
    ok = ~np.isnan(rho)
    r = rho[ok]
    frac_neg = float(np.mean(r < 0))
    med = float(np.median(r))
    bt = stats.binomtest(int((r < 0).sum()), len(r), 0.5, alternative="two-sided")

    rng = np.random.default_rng(seed)
    n_s = Rp.shape[1]
    null = []
    for _ in range(n_perm):
        r2 = row_corr(Rp[:, rng.permutation(n_s)], Re)
        null.append(np.mean(r2[~np.isnan(r2)] < 0))
    null = np.array(null)
    p_perm = (np.sum(np.abs(null - 0.5) >= abs(frac_neg - 0.5)) + 1) / (n_perm + 1)

    print(f"  {label:28s} n={len(r):6,}  负相关 {frac_neg:6.1%}  "
          f"中位ρ {med:+.4f}  二项p {bt.pvalue:9.3g}  "
          f"置换零分布 {null.mean():.3f}±{null.std():.3f}  置换p {p_perm:.4g}")
    return {"type": label, "n": len(r), "frac_neg": frac_neg, "median_rho": med,
            "p_binom": bt.pvalue, "null_mean": null.mean(), "null_sd": null.std(),
            "p_perm": p_perm}


def main():
    print("=" * 78)
    print("G1 阴性对照：非内含子事件是否也呈负相关？（判别生物学 vs 技术耦合）")
    print("=" * 78)

    psi = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t",
                      index_col=0, compression="gzip", low_memory=False)
    meta = psi.iloc[:, :14].copy()
    vals = psi.iloc[:, 14:].copy()
    vals.columns = [barcode(c) for c in vals.columns]
    tumor = [c for c in vals.columns if c[13:15] == "01"]
    vals = vals[tumor].loc[:, ~vals[tumor].columns.duplicated()]
    P = vals.apply(pd.to_numeric, errors="coerce")

    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0, compression="gzip")
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[~expr.index.duplicated()]
    common = [c for c in P.columns if c in set(expr.columns)]
    P, E = P[common], expr[common]
    print(f"样本 {len(common)} 例｜事件 {P.shape[0]:,}\n")

    detect = P.notna().sum(axis=1).values / P.shape[1]
    Psd = P.std(axis=1, skipna=True).values
    gene = meta["GeneSymbol"].astype(str).values
    st = meta["SpliceType"].astype(str).values
    in_expr = np.array([g in set(E.index) for g in gene])
    base = (detect >= MIN_DETECT) & in_expr & (Psd > 0)

    rows = []
    for t in ["IR", "ES", "A3", "A5", "ME"]:
        m = base & (st == t)
        if m.sum() < 200:
            print(f"  {t:28s} n={int(m.sum()):,}  （样本量不足，跳过）")
            continue
        Pk = P.loc[m, :]
        Ek = E.loc[pd.Index(gene[m]), :]
        rows.append(analyze(Pk, Ek, f"{t}（全部）", n_perm=1000))

    # 配对设计：同一基因内 IR 事件 vs ES 事件（控制基因层面混杂）
    print("\n  --- 配对设计：同一基因同时含 IR 与 ES 事件 ---")
    mir = base & (st == "IR")
    mes = base & (st == "ES")
    g_ir = pd.Series(gene[mir], index=P.index[mir])
    g_es = pd.Series(gene[mes], index=P.index[mes])
    shared = sorted(set(g_ir) & set(g_es))
    print(f"      共有基因 {len(shared):,} 个")
    if len(shared) >= 100:
        Pir = P.loc[g_ir[g_ir.isin(shared)].index, :]
        Pes = P.loc[g_es[g_es.isin(shared)].index, :]
        for lab, Px, gx in [("IR", Pir, g_ir), ("ES", Pes, g_es)]:
            Ex = E.loc[pd.Index(gx[Px.index]), :]
            rows.append(analyze(Px, Ex, f"配对基因内·{lab}", n_perm=1000, seed=SEED + 3))

    # IR 内部：按平均 PSI 分层（reads 转移量代理）
    print("\n  --- IR 内部分层：按事件平均 PSI（内含子 reads 占比代理）---")
    m = base & (st == "IR")
    Pk = P.loc[m, :]
    Ek = E.loc[pd.Index(gene[m]), :]
    mpsi = Pk.mean(axis=1).values
    qs = np.quantile(mpsi, [1 / 3, 2 / 3])
    for lab, sel in [("低 PSI 三分之一", mpsi <= qs[0]),
                     ("中 PSI 三分之一", (mpsi > qs[0]) & (mpsi <= qs[1])),
                     ("高 PSI 三分之一", mpsi > qs[1])]:
        s = np.asarray(sel)
        rows.append(analyze(Pk[s], Ek[s], f"IR·{lab}", n_perm=1000, seed=SEED + 7))

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/G1_negative_control_eventtypes.tsv", sep="\t", index=False)
    print(f"\n已输出 {OUT}/G1_negative_control_eventtypes.tsv")


if __name__ == "__main__":
    main()
