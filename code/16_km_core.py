#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16 — 用 Kahles 2018 IR PSI 跑核心分析（burden + 连续轴诊断）。

适配：新数据只有 IR 事件（5 列 meta：event_id/event_type/event_chr/
event_coordinates/alt_region_coordinates），无 SpliceType、无基因注释。

用法:
  python 16_km_core.py GBM
  python 16_km_core.py LGG
输出到 ri_atlas/output_km/{cancer}/ ：
  burden.tsv, pc1.txt, summary.json
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
MIN_DETECT = 0.90
MIN_STD = 0.10


def main():
    cancer = sys.argv[1] if len(sys.argv) > 1 else "GBM"
    psi_path = os.path.join(ROOT, "ri_atlas", "data", f"KM_IR_{cancer}.psi.gz")
    outdir = os.path.join(OUT, cancer)
    os.makedirs(outdir, exist_ok=True)

    # 读 PSI（事件 × 样本，5 列 meta）
    psi = pd.read_csv(psi_path, sep="\t", index_col=0, compression="gzip")
    meta = psi.iloc[:, :4]          # event_type, event_chr, event_coordinates, alt_region_coordinates
    vals = psi.iloc[:, 4:]
    vals = vals.apply(pd.to_numeric, errors="coerce")
    # 样本去重：pandas 读入会对重复列名加 .1 后缀，先还原再按原名去重
    vals.columns = [str(c).split(".")[0] for c in vals.columns]
    vals = vals.loc[:, ~vals.columns.duplicated()]
    print(f"[1] {cancer} 矩阵: {vals.shape[0]} 事件 × {vals.shape[1]} 样本")

    # 检出率过滤（在真实样本集上重算）
    n_valid = vals.notna().sum(axis=1).values
    detect = n_valid / vals.shape[1]
    Psd = vals.std(axis=1, skipna=True).values
    ri_mask = detect >= MIN_DETECT
    hi_mask = ri_mask & (Psd >= MIN_STD)
    print(f"[2] 检出率≥{MIN_DETECT}: {int(ri_mask.sum())}；其中高变异(std≥{MIN_STD}): {int(hi_mask.sum())}")

    # burden = 高检出 IR 事件平均 PSI
    P_ri = vals.loc[ri_mask, :]
    burden = P_ri.mean(axis=0, skipna=True)
    burden.name = "RI_burden"
    print(f"[3] RI burden: 均值 {burden.mean():.4f}，范围 {burden.min():.4f}–{burden.max():.4f}")

    # 连续轴诊断：高变异 IR 事件 PCA，看 PC1 方差占比 + PC1 vs burden 相关
    X = vals.loc[hi_mask, :].T          # 样本 × 高变异事件
    Xz = ((X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)).fillna(0.0)
    pca = PCA(n_components=min(3, Xz.shape[0]))
    pc = pca.fit_transform(Xz.values)
    evr = pca.explained_variance_ratio_
    pc1 = pc[:, 0]
    rho = float(np.corrcoef(pc1, burden.values)[0, 1])
    print(f"[4] PC1 方差占比 {evr[0]*100:.1f}%，PC2 {evr[1]*100:.1f}%")
    print(f"    PC1 vs RI burden 相关 ρ = {rho:.3f}")

    # 保存
    burden.to_frame().to_csv(os.path.join(outdir, "burden.tsv"), sep="\t")
    np.save(os.path.join(outdir, "pc1.npy"), pc1)
    summary = {
        "cancer": cancer,
        "n_samples": int(vals.shape[1]),
        "n_events": int(vals.shape[0]),
        "n_detect": int(ri_mask.sum()),
        "n_hivar": int(hi_mask.sum()),
        "burden_mean": float(burden.mean()),
        "burden_range": [float(burden.min()), float(burden.max())],
        "pc1_var_ratio": float(evr[0]),
        "pc2_var_ratio": float(evr[1]),
        "pc1_vs_burden_rho": rho,
    }
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
    print("[5] 完成 ->", outdir)


if __name__ == "__main__":
    main()
