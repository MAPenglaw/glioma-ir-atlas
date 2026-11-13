#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""32 — 泛癌聚类诊断：检验「IR 是连续轴、无离散亚型」是否为泛癌普遍现象。

对 32 个癌种的每一个单独复现 B 线的诊断流程：
  burden      = 该癌种内检出率 ≥90% 的事件上的平均 PSI（与 B 线 16 号口径一致）
  高变异事件  = 该癌种内检出率 ≥90% 且 PSI 标准差 ≥0.10（与 B 线 22 号口径一致）
  诊断        = PC1 方差解释率 / PC1-burden Spearman / silhouette(k=2) / PAC / GMM BIC

与 22 号的一个差异：22 号在原始特征空间做 k-means 与 GMM；泛癌下特征维度可达数千，
GMM 的 full covariance 必然奇异，因此这里先把样本投到前 min(10, n/4) 个主成分上再诊断。
GBM / LGG 会同时给出结果，可与 22 号的原始空间结果互相印证。

输出（output_pancan/）：
  pancan_diagnostics.tsv   每癌种一行
  pancan_diagnostics.json  同上（含 ΔAUC）
"""
import os
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = os.path.join(ROOT, "output_pancan")

MIN_DETECT = 0.90
MIN_STD = 0.10
MIN_N = 40          # 样本量门槛
N_ITER = 60
SUB_SAMPLE = 0.8
SUB_FEATURE = 0.8
SEED = 20260928

t0 = time.time()

mat = np.load(os.path.join(OUT, "mat_hq.npy"))
samples = pd.read_csv(os.path.join(OUT, "mat_hq_samples.tsv"), sep="\t")
print(f"[load] 矩阵 {mat.shape} | 样本 {len(samples)}", flush=True)
print(f"[nan ] 总缺失率 {np.isnan(mat).mean():.4f}", flush=True)


def diagnostics(S, seed):
    """S: 样本 × PC 得分矩阵。返回 silhouette / PAC / ΔAUC / BIC。"""
    n_s, n_f = S.shape
    rng = np.random.default_rng(seed)
    ks = [2, 3, 4]
    conn = {k: np.zeros((n_s, n_s)) for k in ks}
    for it in range(N_ITER):
        rs = rng.choice(n_s, max(2, int(n_s * SUB_SAMPLE)), replace=False)
        fs = rng.choice(n_f, max(1, int(n_f * SUB_FEATURE)), replace=False)
        sub = S[np.ix_(rs, fs)]
        for k in ks:
            lab_full = np.full(n_s, -1, dtype=int)
            km = KMeans(n_clusters=k, n_init=10, random_state=seed + it)
            lab_full[rs] = km.fit_predict(sub)
            idx = np.where(lab_full >= 0)[0]
            same = (lab_full[idx][:, None] == lab_full[idx][None, :])
            conn[k][np.ix_(idx, idx)] += same

    pac, auc = {}, {}
    grid = np.linspace(0, 1, 101)
    for k in ks:
        C = conn[k] / N_ITER
        np.fill_diagonal(C, 1.0)
        iu = np.triu_indices(n_s, k=1)
        v = C[iu]
        pac[k] = float(np.mean((v > 0.1) & (v < 0.9)))
        vs = np.sort(v)
        cdf = np.searchsorted(vs, grid, side="right") / len(vs)
        auc[k] = float(np.trapezoid(cdf, grid))
    # ΔAUC(k) = AUC(k) - AUC(k-1)；k=1 时共识恒为 1，故 AUC(1) = 1.0
    auc[1] = 1.0
    dauc = {k: float(auc[k] - auc[k - 1]) for k in ks if k - 1 in auc}

    km2 = KMeans(n_clusters=2, n_init=50, random_state=seed)
    lab = km2.fit_predict(S)
    sil = float(silhouette_score(S, lab, sample_size=min(n_s, 1000),
                                 random_state=seed))
    bic = {}
    for k in (1, 2):
        g = GaussianMixture(n_components=k, covariance_type="full",
                            random_state=seed, max_iter=300)
        g.fit(S)
        bic[k] = float(g.bic(S))
    return sil, pac, dauc, bic


rows = []
for proj, grp in samples.groupby("projectID"):
    idx = grp.index.to_numpy()
    n_s = len(idx)
    if n_s < MIN_N:
        print(f"[skip] {proj} n={n_s} < {MIN_N}", flush=True)
        continue
    M = mat[:, idx]                                   # (n_event, n_s)
    det = np.count_nonzero(~np.isnan(M), axis=1) / n_s

    m_b = det >= MIN_DETECT
    burden = np.nanmean(M[m_b], axis=0)

    with np.errstate(all="ignore"):
        sd_ev = np.nanstd(M, axis=1)
    m_v = (det >= MIN_DETECT) & (sd_ev >= MIN_STD)
    n_v = int(m_v.sum())
    if n_v < 20:
        print(f"[skip] {proj} 高变异事件仅 {n_v}", flush=True)
        continue

    Mv = M[m_v]
    mu = np.nanmean(Mv, axis=1, keepdims=True)
    sdev = np.nanstd(Mv, axis=1, keepdims=True)
    X = np.nan_to_num((Mv - mu) / (sdev + 1e-12)).T    # (n_s, n_v)

    n_comp = int(min(10, n_s // 4, X.shape[1]))
    n_comp = max(n_comp, 2)
    pca = PCA(n_components=n_comp, svd_solver="randomized", random_state=SEED)
    S = pca.fit_transform(X)
    pc1_evr = float(pca.explained_variance_ratio_[0])
    rho = float(spearmanr(S[:, 0], burden).statistic)

    sil, pac, dauc, bic = diagnostics(S, SEED)

    rows.append({
        "projectID": proj, "n_samples": n_s,
        "n_burden_events": int(m_b.sum()), "n_hivar_events": n_v,
        "PC1_evr": round(pc1_evr, 4),
        "rho_PC1_burden": round(rho, 4),
        "silhouette_k2": round(sil, 4),
        "PAC_k2": round(pac[2], 4),
        "PAC_k3": round(pac[3], 4),
        "delta_AUC_k2": round(dauc[2], 4),
        "BIC_k1": round(bic[1], 1), "BIC_k2": round(bic[2], 1),
        "BIC_prefers_k1": bic[1] < bic[2],
        "burden_median": round(float(np.median(burden)), 4),
    })
    print(f"[ok] {proj:12s} n={n_s:4d} 变异={n_v:5d} PC1={pc1_evr:.3f} "
          f"rho={rho:.3f} sil={sil:.3f} PAC={pac[2]:.3f} "
          f"BIC1<2={bic[1] < bic[2]}", flush=True)

df = pd.DataFrame(rows).sort_values("n_samples", ascending=False)
df.to_csv(os.path.join(OUT, "pancan_diagnostics.tsv"), sep="\t", index=False)
df.to_json(os.path.join(OUT, "pancan_diagnostics.json"),
           orient="records", indent=1)

print(f"\n=== 汇总（{len(df)} 个癌种）===", flush=True)
print(f"silhouette k=2  中位数 {df.silhouette_k2.median():.3f} "
      f"（范围 {df.silhouette_k2.min():.3f}–{df.silhouette_k2.max():.3f}）", flush=True)
print(f"PAC k=2         中位数 {df.PAC_k2.median():.3f} "
      f"（范围 {df.PAC_k2.min():.3f}–{df.PAC_k2.max():.3f}）", flush=True)
print(f"BIC 偏好 k=1    {int(df.BIC_prefers_k1.sum())} / {len(df)} 个癌种", flush=True)
print(f"PC1 方差解释率  中位数 {df.PC1_evr.median():.3f}", flush=True)
print(f"rho(PC1,burden) 中位数 {df.rho_PC1_burden.median():.3f} "
      f"（最小 {df.rho_PC1_burden.min():.3f}）", flush=True)
print(f"\n[done] 用时 {(time.time()-t0)/60:.1f} min", flush=True)
