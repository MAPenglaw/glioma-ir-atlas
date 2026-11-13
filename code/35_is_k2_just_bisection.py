#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""35 — A-3 补强：k=2 的「亚型」到底是不是把连续轴切成两半？

32 号里 31 个癌种有 16 个的 GMM BIC 偏好 k=2。这本身不能直接读成「存在两个亚型」：
一个连续分布同样可以被两个高斯分量拟合得更好。因此这里直接检验
k=2 的划分与 IR burden 的关系：

  若 k-means(k=2) 的标签几乎完全由 burden 的高低决定（AUC → 1），
  说明所谓的两类只是连续轴的二分，而非离散亚型。

同时报告：把样本按 burden 中位数切成高低两组，与 k=2 标签的一致率（ARI）。

输出：output_pancan/k2_bisection.tsv
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = os.path.join(ROOT, "output_pancan")

MIN_DETECT = 0.90
MIN_STD = 0.10
MIN_N = 40
SEED = 20260928

mat = np.load(os.path.join(OUT, "mat_hq.npy"))
samples = pd.read_csv(os.path.join(OUT, "mat_hq_samples.tsv"), sep="\t")
print(f"[load] {mat.shape}", flush=True)

rows = []
for proj, grp in samples.groupby("projectID"):
    idx = grp.index.to_numpy()
    n_s = len(idx)
    if n_s < MIN_N:
        continue
    M = mat[:, idx]
    det = np.count_nonzero(~np.isnan(M), axis=1) / n_s

    m_b = det >= MIN_DETECT
    burden = np.nanmean(M[m_b], axis=0)

    with np.errstate(all="ignore"):
        sdev = np.nanstd(M, axis=1)
    m_v = (det >= MIN_DETECT) & (sdev >= MIN_STD)
    if m_v.sum() < 20:
        continue
    Mv = M[m_v]
    mu = np.nanmean(Mv, axis=1, keepdims=True)
    sd2 = np.nanstd(Mv, axis=1, keepdims=True)
    X = np.nan_to_num((Mv - mu) / (sd2 + 1e-12)).T

    n_comp = max(2, int(min(10, n_s // 4, X.shape[1])))
    S = PCA(n_components=n_comp, svd_solver="randomized",
            random_state=SEED).fit_transform(X)

    lab = KMeans(n_clusters=2, n_init=50, random_state=SEED).fit_predict(S)

    # k=2 标签区分 burden 的能力（AUC）
    a = burden[lab == 0]
    b = burden[lab == 1]
    U = mannwhitneyu(a, b, alternative="two-sided")
    auc = float(U.statistic / (len(a) * len(b)))
    auc = max(auc, 1 - auc)          # 标签顺序任意，取靠近 1 的一侧

    # 与「burden 中位数二分」的一致率
    med_split = (burden > np.median(burden)).astype(int)
    ari = float(adjusted_rand_score(lab, med_split))

    rows.append({
        "projectID": proj, "n_samples": n_s,
        "AUC_k2label_vs_burden": round(auc, 4),
        "ARI_k2label_vs_median_split": round(ari, 4),
        "k2_cluster_sizes": f"{int((lab==0).sum())}/{int((lab==1).sum())}",
        "balancing": round(min((lab==0).sum(), (lab==1).sum()) / n_s, 3),
        "p_mw": float(U.pvalue),
    })
    print(f"[ok] {proj:12s} n={n_s:4d} AUC={auc:.3f} ARI={ari:.3f} "
          f"簇={int((lab==0).sum())}/{int((lab==1).sum())}", flush=True)

df = pd.DataFrame(rows).sort_values("AUC_k2label_vs_burden", ascending=False)
df.to_csv(os.path.join(OUT, "k2_bisection.tsv"), sep="\t", index=False)

print(f"\n=== 汇总（{len(df)} 个癌种）===", flush=True)
print(f"AUC(k=2 标签 vs burden) 中位数 {df.AUC_k2label_vs_burden.median():.3f} "
      f"（范围 {df.AUC_k2label_vs_burden.min():.3f}–"
      f"{df.AUC_k2label_vs_burden.max():.3f}）", flush=True)
print(f"ARI(k=2 标签 vs 中位数二分) 中位数 "
      f"{df.ARI_k2label_vs_median_split.median():.3f} "
      f"（范围 {df.ARI_k2label_vs_median_split.min():.3f}–"
      f"{df.ARI_k2label_vs_median_split.max():.3f}）", flush=True)
print(f"AUC ≥ 0.9 的癌种: {int((df.AUC_k2label_vs_burden >= 0.9).sum())}"
      f" / {len(df)}", flush=True)
print("\nGBM / LGG 单独看：", flush=True)
print(df[df.projectID.isin(["TCGA-GBM", "TCGA-LGG"])].to_string(index=False),
      flush=True)
print("\n[done]", flush=True)
