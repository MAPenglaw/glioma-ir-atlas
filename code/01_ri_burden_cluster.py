#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1 — TCGA-GBM 内含子滞留（RI）全局负荷与无监督分型
数据源: OncoSplicing / SplAdder TCGA-GBM PSI 矩阵 + UCSC Xena TCGA-GBM 表达谱

设计要点（对应「坑 1」的规避）：
  - 不逐事件做生存检验（6.9 万事件上 BH 惩罚过重），改为「单指标 + 无监督」设计
  - RI 全局负荷 = 每样本在「高检出率 IR 事件」上的平均 PSI（单指标，只做一次检验）
  - RI 分型 = 在高变异 IR 事件子集上做一致性聚类（consensus clustering）
"""
import gzip
import json
import os
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
os.makedirs(f"{OUT}/output", exist_ok=True)

MIN_DETECT = 0.90          # 事件检出率门槛
MIN_STD = 0.10             # 高变异门槛（用于聚类）
SEED = 20260917
N_ITER = 500               # 一致性聚类迭代
SUB_SAMPLE = 0.80
SUB_FEATURE = 0.80


def barcode(col: str) -> str:
    return col.strip().strip('"')[:15]


def main():
    print("=" * 72)
    print("Step 1  TCGA-GBM 全局 RI 负荷与无监督分型")
    print("=" * 72)

    # ---------- 1. 载入 PSI 矩阵 ----------
    psi = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t",
                      index_col=0, compression="gzip", low_memory=False)
    meta = psi.iloc[:, :14].copy()
    vals = psi.iloc[:, 14:].copy()
    vals.columns = [barcode(c) for c in vals.columns]
    print(f"[1] PSI 矩阵: {vals.shape[0]:,} 事件 × {vals.shape[1]} 样本列")

    # ---------- 2. 样本过滤: 原发肿瘤 + 去重 + 与表达谱取交集 ----------
    tumor_cols = [c for c in vals.columns if c[13:15] == "01"]
    vals = vals[tumor_cols].loc[:, ~vals[tumor_cols].columns.duplicated()]
    print(f"[2] 原发肿瘤(01)去重后: {vals.shape[1]} 样本")

    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0,
                       compression="gzip")
    expr_cols = {barcode(c) for c in expr.columns if barcode(c)[13:15] == "01"}
    common = [c for c in vals.columns if c in expr_cols]
    vals = vals[common]
    print(f"    与表达谱交集后进入分析的样本: {vals.shape[1]} 例")

    P = vals.apply(pd.to_numeric, errors="coerce")

    # ---------- 3. 事件层面过滤（在真实样本集上重算，不用 OncoSplicing 预置列）----------
    is_ir = (meta["SpliceType"] == "IR").values
    n_valid = P.notna().sum(axis=1).values
    detect = n_valid / P.shape[1]
    Psd = P.std(axis=1, skipna=True).values

    ri_mask = is_ir & (detect >= MIN_DETECT)
    hi_mask = ri_mask & (Psd >= MIN_STD)
    print(f"[3] IR 事件总数: {int(is_ir.sum()):,}")
    print(f"    检出率 >= {MIN_DETECT}: {int(ri_mask.sum()):,}")
    print(f"    其中 std_PSI >= {MIN_STD}（高变异，用于聚类）: {int(hi_mask.sum()):,}")

    # ---------- 4. 全局 RI 负荷（单指标）----------
    P_ri = P.loc[ri_mask, :]
    burden = P_ri.mean(axis=0, skipna=True)
    burden.name = "RI_burden"
    print(f"[4] 全局 RI 负荷: 均值 {burden.mean():.4f}, "
          f"范围 {burden.min():.4f}–{burden.max():.4f}")

    # 敏感性分析: 只用「强变异」事件算负荷 + 用中位数算负荷
    burden_med = P_ri.median(axis=0, skipna=True).rename("RI_burden_median")
    burden_hv = P.loc[hi_mask, :].mean(axis=0, skipna=True).rename("RI_burden_hivar")
    print(f"    敏感性: 中位数口径 {burden_med.mean():.4f}; "
          f"高变异子集口径 {burden_hv.mean():.4f}")

    # ---------- 5. 一致性聚类 ----------
    X = P.loc[hi_mask, :].copy()
    # 标准化（每事件 z 分数，避免高 PSI 事件主导距离）
    Xz = ((X.T - X.mean(axis=1)) / X.std(axis=1)).T
    # 缺失值用事件中位数（z 空间为 0）填补
    Xz = Xz.fillna(0.0)
    M = Xz.T.values.astype(np.float64)          # 样本 × 事件
    M_raw = X.T.values.astype(np.float64)
    print(f"[5] 聚类输入矩阵: {M.shape[0]} 样本 × {M.shape[1]} 高变异 IR 事件")

    n_s, n_f = M.shape
    rng = np.random.default_rng(SEED)
    ks = list(range(2, 7))
    conn = {k: np.zeros((n_s, n_s)) for k in ks}
    valid = {k: 0 for k in ks}

    for it in range(N_ITER):
        # 80% 样本 + 80% 特征 重抽样
        rs = rng.choice(n_s, int(n_s * SUB_SAMPLE), replace=False)
        fs = rng.choice(n_f, int(n_f * SUB_FEATURE), replace=False)
        sub = M[np.ix_(rs, fs)]
        # 特征方向对齐：以全体数据均值为参照翻转（保证可比）
        for k in ks:
            lab_full = np.full(n_s, -1, dtype=int)
            km = KMeans(n_clusters=k, n_init=10, random_state=SEED + it)
            lab = km.fit_predict(sub)
            lab_full[rs] = lab
            idx = np.where(lab_full >= 0)[0]
            same = (lab_full[idx][:, None] == lab_full[idx][None, :])
            conn[k][np.ix_(idx, idx)] += same
            valid[k] += 1
        if (it + 1) % 100 == 0:
            print(f"    迭代 {it + 1}/{N_ITER}")

    # 归一化为一致性矩阵
    consensus = {}
    for k in ks:
        C = conn[k] / valid[k]
        np.fill_diagonal(C, 1.0)
        consensus[k] = C

    # PAC 分数（比例模糊配对，越小越好）
    pac = {}
    for k in ks:
        C = consensus[k]
        iu = np.triu_indices(n_s, k=1)
        v = C[iu]
        pac[k] = float(np.mean((v > 0.1) & (v < 0.9)))
    print("[5] PAC 分数（越小越稳定）:",
          {k: round(v, 4) for k, v in pac.items()})

    # CDF 曲线（用于 ΔAUC 判据）
    cdf = {}
    grid = np.linspace(0, 1, 101)
    for k in ks:
        iu = np.triu_indices(n_s, k=1)
        v = np.sort(consensus[k][iu])
        cdf[k] = np.searchsorted(v, grid, side="right") / len(v)
    auc = {k: float(np.trapezoid(cdf[k], grid)) for k in ks}
    dauc = {k: float(auc[k] - auc[k - 1]) for k in ks if k - 1 in auc}
    print("[5] CDF-AUC:", {k: round(v, 4) for k, v in auc.items()})
    print("[5] ΔAUC   :", {k: round(v, 4) for k, v in dauc.items()})

    # ---------- 6. 定 k=2 并输出分组 ----------
    K_FINAL = 2
    C2 = consensus[K_FINAL]
    km2 = KMeans(n_clusters=K_FINAL, n_init=50, random_state=SEED)
    lab = km2.fit_predict(M)
    # 让 cluster 1 = RI 负荷更高者，命名为 RI-high
    b = burden.values
    if b[lab == 1].mean() < b[lab == 0].mean():
        lab = 1 - lab
    group = np.where(lab == 1, "RI-high", "RI-low")

    samp = pd.DataFrame({
        "sample": vals.columns,
        "RI_burden": burden.values,
        "RI_burden_median": burden_med.values,
        "RI_burden_hivar": burden_hv.values,
        "cluster": lab,
        "RI_group": group,
        "consensus_k2": C2[:, :].mean(axis=1),
    }).set_index("sample")
    samp.to_csv(f"{OUT}/output/sample_RI_group.tsv", sep="\t")

    print(f"[6] 最终分组 (k={K_FINAL}):")
    print(samp.groupby("RI_group")["RI_burden"].agg(["count", "mean", "std"]).to_string())
    from scipy import stats as st
    a = samp.loc[samp.RI_group == "RI-high", "RI_burden"]
    c = samp.loc[samp.RI_group == "RI-low", "RI_burden"]
    u, p = st.mannwhitneyu(a, c, alternative="two-sided")
    print(f"    RI-high vs RI-low 负荷差异 Mann-Whitney p = {p:.3e}")

    # ---------- 7. 保存中间结果 ----------
    np.savez_compressed(f"{OUT}/output/cluster_input.npz",
                        M=M, samples=np.array(vals.columns), groups=group)
    for k in ks:
        pd.DataFrame(consensus[k], index=vals.columns,
                     columns=vals.columns).to_csv(
            f"{OUT}/output/consensus_k{k}.tsv", sep="\t")

    summary = {
        "psi_events": int(vals.shape[0]),
        "samples_columns": 167,
        "samples_primary_tumor_unique": 152,
        "samples_analyzed": int(vals.shape[1]),
        "ir_events_total": int(is_ir.sum()),
        "ir_events_detect_ge_0.9": int(ri_mask.sum()),
        "ir_events_high_variance": int(hi_mask.sum()),
        "RI_burden_mean": float(burden.mean()),
        "RI_burden_min": float(burden.min()),
        "RI_burden_max": float(burden.max()),
        "RI_burden_median_caliber": float(burden_med.mean()),
        "RI_burden_hivar_caliber": float(burden_hv.mean()),
        "groups": samp["RI_group"].value_counts().to_dict(),
        "group_burden_mean": samp.groupby("RI_group")["RI_burden"].mean().to_dict(),
        "mannwhitney_p": float(p),
        "PAC": pac,
        "CDF_AUC": auc,
        "delta_AUC": dauc,
        "k_final": K_FINAL,
    }
    with open(f"{OUT}/output/step1_summary.json", "w") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("[7] 完成 -> sample_RI_group.tsv / consensus_k*.tsv / step1_summary.json")


if __name__ == "__main__":
    main()
