#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""22 — 用 Kahles 2018 数据重算聚类诊断（silhouette / PAC / ΔAUC / GMM BIC）。

旧版 01 脚本在 OncoSplicing 高变异事件上算的诊断值（silhouette 0.143 / PAC 0.525 等）
已随数据源切换失效。这里用 Kahles 的 KM_IR_{GBM,LGG}.psi.gz 高变异事件重算，
方法与旧版一致：consensus clustering（80% 样本 × 80% 特征重采样）。
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
DATA = os.path.join(ROOT, "ri_atlas", "data")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")

MIN_DETECT = 0.90
MIN_STD = 0.10
N_ITER = 100
SUB_SAMPLE = 0.8
SUB_FEATURE = 0.8
SEED = 20260927


def load_hivar(cancer):
    psi = pd.read_csv(os.path.join(DATA, f"KM_IR_{cancer}.psi.gz"), sep="\t",
                      compression="gzip", low_memory=False)
    vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    vals = vals.loc[:, ~vals.columns.duplicated()]
    detect = vals.notna().mean(axis=1)
    std = vals.std(axis=1, skipna=True)
    mask = (detect >= MIN_DETECT) & (std >= MIN_STD)
    M = vals.loc[mask]
    Mv = M.values.astype(float)
    mu = np.nanmean(Mv, axis=1, keepdims=True)
    sd = np.nanstd(Mv, axis=1, keepdims=True)
    Xz = (Mv - mu) / (sd + 1e-12)
    X = np.nan_to_num(Xz, nan=0.0).T
    return X, int(mask.sum())


def diagnostics(X, n_s, n_f):
    rng = np.random.default_rng(SEED)
    ks = [2, 3, 4, 5, 6]
    conn = {k: np.zeros((n_s, n_s)) for k in ks}
    valid = {k: 0 for k in ks}
    for it in range(N_ITER):
        rs = rng.choice(n_s, int(n_s * SUB_SAMPLE), replace=False)
        fs = rng.choice(n_f, int(n_f * SUB_FEATURE), replace=False)
        sub = X[np.ix_(rs, fs)]
        for k in ks:
            lab_full = np.full(n_s, -1, dtype=int)
            km = KMeans(n_clusters=k, n_init=10, random_state=SEED + it)
            lab = km.fit_predict(sub)
            lab_full[rs] = lab
            idx = np.where(lab_full >= 0)[0]
            same = (lab_full[idx][:, None] == lab_full[idx][None, :])
            conn[k][np.ix_(idx, idx)] += same
            valid[k] += 1

    consensus = {}
    for k in ks:
        C = conn[k] / max(valid[k], 1)
        np.fill_diagonal(C, 1.0)
        consensus[k] = C

    # PAC
    pac = {}
    for k in ks:
        iu = np.triu_indices(n_s, k=1)
        v = consensus[k][iu]
        pac[k] = float(np.mean((v > 0.1) & (v < 0.9)))

    # ΔAUC
    cdf, auc = {}, {}
    grid = np.linspace(0, 1, 101)
    for k in ks:
        iu = np.triu_indices(n_s, k=1)
        v = np.sort(consensus[k][iu])
        cdf[k] = np.searchsorted(v, grid, side="right") / len(v)
        auc[k] = float(np.trapezoid(cdf[k], grid))
    dauc = {k: float(auc[k] - auc[k - 1]) for k in ks if k - 1 in auc}

    # silhouette k=2
    km2 = KMeans(n_clusters=2, n_init=50, random_state=SEED)
    lab = km2.fit_predict(X)
    sil = float(silhouette_score(X, lab, sample_size=min(n_s, 1000)))

    # GMM BIC k=1 vs k=2
    bic = {}
    for k in (1, 2):
        g = GaussianMixture(n_components=k, covariance_type="full", random_state=SEED)
        g.fit(X)
        bic[k] = float(g.bic(X))
    return {"silhouette_k2": sil, "PAC": {str(k): pac[k] for k in ks},
            "delta_AUC": {str(k): dauc.get(k) for k in ks},
            "BIC": {str(k): bic[k] for k in bic},
            "BIC_k1_k2": f"{bic[1]:.1f} vs {bic[2]:.1f} (lower=better)"}


if __name__ == "__main__":
    result = {}
    for cancer in ("GBM", "LGG"):
        X, n_hivar = load_hivar(cancer)
        n_s, n_f = X.shape
        print(f"[{cancer}] 高变异事件 {n_hivar}，样本 {n_s} × 特征 {n_f}", flush=True)
        result[cancer] = diagnostics(X, n_s, n_f)
        print(f"  silhouette k=2: {result[cancer]['silhouette_k2']:.4f}")
        print(f"  PAC: {result[cancer]['PAC']}")
        print(f"  ΔAUC: {result[cancer]['delta_AUC']}")
        print(f"  BIC: {result[cancer]['BIC_k1_k2']}")
    json.dump(result, open(os.path.join(OUT, "cluster_diagnostics.json"), "w"), indent=1)
    print("DONE -> output_km/cluster_diagnostics.json")
