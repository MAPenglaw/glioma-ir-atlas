#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1c — 判别实验：PSI–mRNA 负耦合是否依赖 NMD 机器？

逻辑：
  · 若耦合来自「内含子 reads 挤占外显子定量」的技术效应 → 与 NMD 机器表达无关
  · 若耦合来自 NMD 降解 → NMD 机器表达高的样本，耦合应更强

设计（调节分析，避免分组造成的 range restriction）：
  对每个 IR 事件拟合  rank(mRNA) ~ rank(PSI) + rank(PSI) × NMD_score
  检验交互项 beta_int 的分布是否显著偏向负（NMD 强化负耦合）
  置换零分布：打乱样本的 NMD 评分标签
"""
import gc
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
BASE = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = f"{BASE}/output"
MIN_DETECT = 0.90
N_PERM = 500
SEED = 20260917

NMD_MACHINERY = ["UPF1", "UPF2", "UPF3A", "UPF3B", "SMG1", "SMG5", "SMG6",
                 "SMG7", "SMG8", "SMG9", "EIF4A3", "MAGOH", "RBM8A", "CASC3",
                 "NBAS", "DHX34", "PNRC2", "WIBG", "GSPT1", "ETF1"]


def barcode(c):
    return str(c).strip().strip('"')[:15]


def rank_rows(A):
    out = np.empty_like(A)
    for i in range(A.shape[0]):
        out[i] = rankdata(A[i])
    return out


def interaction_beta(Xr, Yr, z):
    """逐事件 OLS：y ~ x + x*z  （x、y 已秩化并中心化），返回 beta_int"""
    n = Xr.shape[1]
    zc = z - z.mean()
    xc = Xr - Xr.mean(axis=1, keepdims=True)
    yc = Yr - Yr.mean(axis=1, keepdims=True)
    xz = xc * zc
    # 逐行解 3×3 正规方程
    b = np.empty(Xr.shape[0])
    for i in range(Xr.shape[0]):
        A = np.column_stack([xc[i], zc, xz[i]])
        try:
            coef, *_ = np.linalg.lstsq(A, yc[i], rcond=None)
            b[i] = coef[2]
        except np.linalg.LinAlgError:
            b[i] = np.nan
    return b


def main():
    print("=" * 78)
    print("G1c  耦合是否依赖 NMD 机器？（判别生物学 vs 技术耦合）")
    print("=" * 78)

    meta = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t", index_col=0,
                       compression="gzip", usecols=[0, 1, 4], low_memory=False)
    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0, compression="gzip")
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[~expr.index.duplicated()]

    rows = []
    reader = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t", index_col=0,
                         compression="gzip", chunksize=8000, low_memory=False)
    for ch in reader:
        v = ch.iloc[:, 13:]
        v.columns = [barcode(c) for c in v.columns]
        tum = [c for c in v.columns if len(c) >= 15 and c[13:15] == "01"]
        v = v[tum].loc[:, ~v[tum].columns.duplicated()]
        v = v.apply(pd.to_numeric, errors="coerce")
        v = v[[c for c in v.columns if c in set(expr.columns)]]
        rows.append(v.loc[v.index.intersection(meta.index)])
        del ch, v
        gc.collect()
    P = pd.concat(rows)
    del rows
    gc.collect()
    print(f"[1] PSI {P.shape[0]:,} 事件 × {P.shape[1]} 样本")

    detect = P.notna().sum(1).values / P.shape[1]
    Psd = P.std(axis=1, skipna=True).values
    gene = meta.loc[P.index, "GeneSymbol"].astype(str).values
    st = meta.loc[P.index, "SpliceType"].astype(str).values
    base = (detect >= MIN_DETECT) & np.isin(gene, np.asarray(expr.index)) & (Psd > 0)
    m = base & (st == "IR")
    Pk = P.loc[m, :]
    print(f"[2] 可分析 IR 事件 {Pk.shape[0]:,}")

    # NMD 机器评分
    nmd_g = [g for g in NMD_MACHINERY if g in set(expr.index)]
    print(f"[3] NMD 机器可用基因 {len(nmd_g)} 个：{', '.join(nmd_g)}")
    Z = expr.loc[nmd_g, Pk.columns].rank(axis=1)
    z = Z.mean(axis=0)
    z = (z - z.mean()) / z.std()
    print(f"    NMD 评分范围 {z.min():+.2f} – {z.max():+.2f}")

    # 秩化
    med = Pk.median(axis=1).values
    A = Pk.values.copy()
    ip = np.where(np.isnan(A))
    A[ip] = np.take(med, ip[0])
    Xr = rank_rows(A)
    Yr = rank_rows(expr.loc[gene[m], Pk.columns].values)
    print("[4] 拟合交互模型 y ~ x + x*z ...")
    beta = interaction_beta(Xr, Yr, z.values)
    b_ok = beta[~np.isnan(beta)]
    print(f"    交互系数：均值 {b_ok.mean():+.5f}｜中位 {np.median(b_ok):+.5f}"
          f"｜负值占比 {np.mean(b_ok < 0):.1%}")

    # 置换
    print(f"[5] 置换零分布（{N_PERM} 次，打乱 NMD 评分标签）...")
    rng = np.random.default_rng(SEED)
    n = len(z)
    null = []
    for _ in range(N_PERM):
        null.append(np.mean(interaction_beta(Xr, Yr, z.values[rng.permutation(n)])))
    null = np.array(null)
    obs = b_ok.mean()
    p_perm = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
    print(f"    零分布均值 {null.mean():+.5f} ± {null.std():.5f}"
          f"｜观测值 {obs:+.5f}｜置换 p = {p_perm:.4g}")

    # 对照：用完全无关的评分（细胞周期无关基因）做负向对照
    ctrl_g = [g for g in ["HBB", "ALB", "INS", "MYH7", "KRT1"] if g in set(expr.index)]
    if ctrl_g:
        zc = expr.loc[ctrl_g, Pk.columns].rank(axis=1).mean(axis=0)
        zc = (zc - zc.mean()) / zc.std()
        bc = interaction_beta(Xr, Yr, zc.values)
        bc = bc[~np.isnan(bc)]
        print(f"[6] 阴性对照（{','.join(ctrl_g)}）：交互均值 {bc.mean():+.5f}"
              f"（应接近 0）")

    pd.DataFrame({"event": Pk.index, "gene": gene[m], "beta_int": beta}).to_csv(
        f"{OUT}/G1c_nmd_interaction.tsv", sep="\t", index=False)
    pd.DataFrame({"metric": ["n_events", "beta_mean", "beta_median", "frac_neg",
                             "null_mean", "null_sd", "p_perm"],
                  "value": [len(b_ok), obs, float(np.median(b_ok)),
                            float(np.mean(b_ok < 0)), float(null.mean()),
                            float(null.std()), p_perm]}).to_csv(
        f"{OUT}/G1c_summary.tsv", sep="\t", index=False)
    print(f"[7] 已输出 {OUT}/G1c_*.tsv")


if __name__ == "__main__":
    main()
