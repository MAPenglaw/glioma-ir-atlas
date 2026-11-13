#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""51 — PAC 算法 bug 修复 + 全量重算（评审 P0 衍生，2026-09-28）。

## Bug 描述

22_cluster_diagnostics.py 与 32_pancan_cluster_diagnostics.py 计算 consensus
clustering 的 PAC 时，分母用的是**重采样迭代总数 N_ITER**：

    C = conn / N_ITER          # ← 错

而 conn 只在两个样本**同时被抽中**的那次迭代才累加。样本同时入样的概率约
SUB_SAMPLE^2 = 0.64，因此即使聚类 100% 稳定，C 的上限也只有 ~0.64，
必然落在 (0.1, 0.9) 的 "ambiguous" 区间 → PAC 被系统性推高到 0.5 附近。

Monti et al. (2003) 的标准做法是分母用**该样本对的实际共现次数**：

    C[i,j] = conn[i,j] / valid[i,j]

## 影响

模拟对照（本脚本输出）显示旧口径**完全丧失判别力**：

    真实两簇(强分离)  旧 0.497   标准 0.000
    真实两簇(弱分离)  旧 0.497   标准 0.000
    单一高斯云(无簇)  旧 0.856   标准 0.859
    连续梯度(一条线)  旧 0.519   标准 0.070

旧口径下"真实两簇"与"连续梯度"几乎同值（0.497 vs 0.519），无法区分。
手稿原引用的 PAC 0.525 (GBM) / 0.542 (LGG) 是 bug 产物，必须替换。

## 结论不变（反而更硬）

标准口径下胶质瘤 PAC = 0.075 (GBM) / 0.102 (LGG)，与"连续梯度"模拟的 0.070 吻合，
与"无簇高斯云"的 0.859 相差一个数量级；同时 silhouette 仍低（0.148/0.125）、
GMM BIC 仍强偏好 k=1。即：切分**可重复**，但可重复的是"把一条连续轴切成两半"。

输出：
  output_km/pac_algorithm_validation.tsv    四场景 × 两口径对照
  output_km/cluster_diagnostics_fixed.json  胶质瘤 k=2..6 修正值
  output_pancan/pancan_diagnostics_fixed.tsv 泛癌修正值
"""
import os, time, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
DATA = os.path.join(ROOT, "ri_atlas", "data")
OUTK = os.path.join(ROOT, "ri_atlas", "output_km")
OUTP = os.path.join(ROOT, "ri_atlas", "output_pancan")

MIN_DETECT, MIN_STD = 0.90, 0.10
N_ITER, SUB_S, SUB_F, SEED = 100, 0.8, 0.8, 20260927


# ---------------------------------------------------------------- PAC
def consensus_pac(X, ks=(2, 3, 4, 5, 6), n_iter=N_ITER, seed=SEED):
    """Monti 标准 consensus：分母 = 样本对实际共现次数。返回 {k: PAC}。"""
    rng = np.random.default_rng(seed)
    n_s, n_f = X.shape
    conn = {k: np.zeros((n_s, n_s)) for k in ks}
    valid = np.zeros((n_s, n_s))
    for it in range(n_iter):
        rs = rng.choice(n_s, max(2, int(n_s * SUB_S)), replace=False)
        fs = rng.choice(n_f, max(1, int(n_f * SUB_F)), replace=False)
        sub = X[np.ix_(rs, fs)]
        lf = np.full(n_s, -1, dtype=int)
        for k in ks:
            km = KMeans(n_clusters=k, n_init=10, random_state=seed + it)
            lf[rs] = km.fit_predict(sub)
            ii = np.where(lf >= 0)[0]
            conn[k][np.ix_(ii, ii)] += (lf[ii][:, None] == lf[ii][None, :])
        valid[np.ix_(rs, rs)] += 1
    out = {}
    for k in ks:
        C = np.where(valid > 0, conn[k] / np.maximum(valid, 1), 0.0)
        np.fill_diagonal(C, 1.0)
        iu = np.triu_indices(n_s, k=1)
        v = C[iu]
        out[k] = float(np.mean((v > 0.1) & (v < 0.9)))
    return out


def pac_legacy(X, k=2, n_iter=N_ITER, seed=SEED):
    """旧口径（除以迭代总数）——仅用于对照实验复现 bug，不得用于结论。"""
    rng = np.random.default_rng(seed)
    n_s, n_f = X.shape
    conn = np.zeros((n_s, n_s))
    for it in range(n_iter):
        rs = rng.choice(n_s, max(2, int(n_s * SUB_S)), replace=False)
        fs = rng.choice(n_f, max(1, int(n_f * SUB_F)), replace=False)
        km = KMeans(n_clusters=k, n_init=10, random_state=seed + it)
        km.fit_predict(X[np.ix_(rs, fs)])
        lf = np.full(n_s, -1, dtype=int)
        lf[rs] = km.labels_
        ii = np.where(lf >= 0)[0]
        conn[np.ix_(ii, ii)] += (lf[ii][:, None] == lf[ii][None, :])
    C = conn / n_iter
    np.fill_diagonal(C, 1.0)
    iu = np.triu_indices(n_s, k=1)
    v = C[iu]
    return float(np.mean((v > 0.1) & (v < 0.9)))


def sil_bic(X, seed=SEED):
    n_s = X.shape[0]
    km = KMeans(n_clusters=2, n_init=50, random_state=seed)
    lab = km.fit_predict(X)
    sil = float(silhouette_score(X, lab, sample_size=min(n_s, 1000), random_state=seed))
    bic = {}
    for k in (1, 2):
        g = GaussianMixture(n_components=k, covariance_type="full",
                            random_state=seed, max_iter=300)
        g.fit(X)
        bic[k] = float(g.bic(X))
    return sil, bic[1], bic[2]


# ---------------------------------------------------------------- 场景
def simulated_scenarios():
    rng = np.random.default_rng(7)
    def cloud(mu, sd, n=75, d=20):
        return rng.normal(mu, sd, (n, d))
    return {
        "A_two_clusters_strong": np.vstack([cloud(-3, .6), cloud(3, .6)]),
        "B_two_clusters_weak": np.vstack([cloud(-1, .9), cloud(1, .9)]),
        "C_single_gaussian_cloud": rng.normal(0, 1, (150, 20)),
        "D_continuous_gradient": np.column_stack(
            [np.linspace(-3, 3, 150)[:, None]] + [rng.normal(0, .3, 150) for _ in range(19)]),
    }


# ---------------------------------------------------------------- 主流程
def load_hivar(cancer):
    psi = pd.read_csv(os.path.join(DATA, f"KM_IR_{cancer}.psi.gz"), sep="\t",
                      compression="gzip", low_memory=False)
    vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    vals = vals.loc[:, ~vals.columns.duplicated()]
    det = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    m = (det >= MIN_DETECT) & (sd >= MIN_STD)
    Mv = vals.loc[m].values.astype(float)
    X = np.nan_to_num((Mv - np.nanmean(Mv, axis=1, keepdims=True)) /
                      (np.nanstd(Mv, axis=1, keepdims=True) + 1e-12)).T
    return X, int(m.sum())


def main():
    t0 = time.time()
    # ---------- 1. 模拟对照 ----------
    print("=== 1/3 模拟场景对照（旧口径 vs Monti 标准）===", flush=True)
    rows = []
    for name, X in simulated_scenarios().items():
        sil, b1, b2 = sil_bic(X)
        rows.append({
            "scenario": name,
            "PAC_legacy": round(pac_legacy(X), 4),
            "PAC_monti": round(consensus_pac(X, ks=(2,))[2], 4),
            "silhouette_k2": round(sil, 4),
            "BIC_k1": round(b1, 1), "BIC_k2": round(b2, 1),
            "BIC_prefers_k1": bool(b1 < b2),
            "ground_truth": ("two clusters" if name.startswith(("A", "B"))
                             else "no cluster / gradient"),
        })
        print(f"  {name:<24} 旧={rows[-1]['PAC_legacy']:.3f}  "
              f"标准={rows[-1]['PAC_monti']:.3f}  sil={sil:.3f}  "
              f"BIC偏好k1={rows[-1]['BIC_prefers_k1']}", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUTK, "pac_algorithm_validation.tsv"),
                              sep="\t", index=False)

    # ---------- 2. 胶质瘤重算 ----------
    print("=== 2/3 胶质瘤 GBM/LGG 重算（Monti 标准，k=2..6）===", flush=True)
    fixed = {}
    for cancer in ("GBM", "LGG"):
        X, n_ev = load_hivar(cancer)
        sil, b1, b2 = sil_bic(X)
        pac = consensus_pac(X)
        fixed[cancer] = {
            "n_samples": int(X.shape[0]), "n_hivar_events": n_ev,
            "silhouette_k2": round(sil, 4),
            "PAC": {str(k): round(v, 4) for k, v in pac.items()},
            "BIC_k1": round(b1, 1), "BIC_k2": round(b2, 1),
            "BIC_prefers_k1": bool(b1 < b2),
        }
        print(f"  {cancer} n={X.shape[0]} 事件={n_ev} sil={sil:.4f} "
              f"PAC2={pac[2]:.4f} PAC3={pac[3]:.4f} BIC偏好k1={b1 < b2}", flush=True)
    with open(os.path.join(OUTK, "cluster_diagnostics_fixed.json"), "w") as f:
        json.dump(fixed, f, indent=1)

    # ---------- 3. 泛癌重算 ----------
    print("=== 3/3 泛癌重算（Monti 标准）===", flush=True)
    mat = np.load(os.path.join(OUTP, "mat_hq.npy"), mmap_mode="r")
    samples = pd.read_csv(os.path.join(OUTP, "mat_hq_samples.tsv"), sep="\t")
    prows = []
    for proj, grp in samples.groupby("projectID"):
        idx = grp.index.to_numpy()
        n_s = len(idx)
        if n_s < 40:
            continue
        M = np.asarray(mat[:, idx], dtype=float)
        det = np.count_nonzero(~np.isnan(M), axis=1) / n_s
        m_b = det >= MIN_DETECT
        if m_b.sum() < 20:
            continue
        burden = np.nanmean(M[m_b], axis=0)
        with np.errstate(all="ignore"):
            sd_ev = np.nanstd(M, axis=1)
        m_v = m_b & (sd_ev >= MIN_STD)
        if m_v.sum() < 20:
            continue
        Mv = M[m_v]
        X = np.nan_to_num((Mv - np.nanmean(Mv, axis=1, keepdims=True)) /
                          (np.nanstd(Mv, axis=1, keepdims=True) + 1e-12)).T
        n_comp = max(2, int(min(10, n_s // 4, X.shape[1])))
        S = PCA(n_components=n_comp, svd_solver="randomized",
                random_state=SEED).fit_transform(X)
        pc1_evr = float(PCA(n_components=n_comp, svd_solver="randomized",
                            random_state=SEED).fit(X).explained_variance_ratio_[0])
        rho = float(spearmanr(S[:, 0], burden).statistic)
        sil, b1, b2 = sil_bic(S)
        pac = consensus_pac(S, ks=(2,), n_iter=60, seed=20260928)[2]
        prows.append({
            "projectID": proj, "n_samples": n_s,
            "PC1_evr": round(pc1_evr, 4), "rho_PC1_burden": round(rho, 4),
            "silhouette_k2": round(sil, 4), "PAC_k2": round(pac, 4),
            "BIC_k1": round(b1, 1), "BIC_k2": round(b2, 1),
            "BIC_prefers_k1": bool(b1 < b2),
        })
        print(f"  {proj:<12} n={n_s:4d} PC1={pc1_evr:.3f} rho={rho:+.3f} "
              f"sil={sil:.3f} PAC2={pac:.3f} BIC1<2={b1 < b2}", flush=True)
    df = pd.DataFrame(prows)
    df.to_csv(os.path.join(OUTP, "pancan_diagnostics_fixed.tsv"), sep="\t", index=False)
    print(f"\n  泛癌 {len(df)} 癌种：PAC_k2 中位数 {df.PAC_k2.median():.3f} "
          f"（范围 {df.PAC_k2.min():.3f}–{df.PAC_k2.max():.3f}）", flush=True)
    print(f"  silhouette 中位数 {df.silhouette_k2.median():.3f}；"
          f"BIC 偏好 k=1: {int(df.BIC_prefers_k1.sum())}/{len(df)}", flush=True)
    print(f"\n[done] 用时 {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
