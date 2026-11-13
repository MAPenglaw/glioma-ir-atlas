#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2a — RI 分型的结构诊断 + 临床特征合并
回答一个必须先回答的问题：这套 RI 矩阵里到底有没有「离散亚型」？

判据：
  1) 一致性矩阵 PAC 分数 vs 置换零分布（随机数据基线）
  2) PCA 主成分方差解释比例（是否存在单一主导轴）
  3) 轮廓系数 silhouette（k=2..6）
  4) 负荷的分布形态（偏态/多峰性检验）
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
MATRIX = "/Users/peng/workspace/projects/project_002_TCGA_CGGA胶质瘤数据库/data/TCGA/GBM_clinicalMatrix.txt"
SURV = f"{OUT}/data/GBM_survival.txt"
SEED = 20260917


def barcode(c):
    return str(c).strip().strip('"')[:15]


def pac_score(M, k, n_iter=200, rng=None):
    n_s = M.shape[0]
    conn = np.zeros((n_s, n_s))
    for it in range(n_iter):
        rs = rng.choice(n_s, int(n_s * 0.8), replace=False)
        fs = rng.choice(M.shape[1], int(M.shape[1] * 0.8), replace=False)
        lab = KMeans(k, n_init=5, random_state=SEED + it).fit_predict(M[np.ix_(rs, fs)])
        same = (lab[:, None] == lab[None, :])
        conn[np.ix_(rs, rs)] += same
    C = conn / n_iter
    iu = np.triu_indices(n_s, 1)
    v = C[iu]
    return float(np.mean((v > 0.1) & (v < 0.9)))


def main():
    rng = np.random.default_rng(SEED)
    d = np.load(f"{OUT}/output/cluster_input.npz", allow_pickle=True)
    M, samples, groups = d["M"], d["samples"], d["groups"]
    samp = pd.read_csv(f"{OUT}/output/sample_RI_group.tsv", sep="\t", index_col=0)
    n_s = M.shape[0]
    print(f"输入: {n_s} 样本 × {M.shape[1]} 高变异 IR 事件")

    res = {}

    # ---------- 1. PCA 结构 ----------
    Z = (M - M.mean(0)) / M.std(0)
    pca = PCA(n_components=20, random_state=SEED)
    pcs = pca.fit_transform(Z)
    evr = pca.explained_variance_ratio_
    print("\n[1] PCA 方差解释（前 10 个主成分）:")
    print("    " + ", ".join(f"PC{i+1}={v*100:.1f}%" for i, v in enumerate(evr[:10])))
    res["pca_evr_top10"] = [float(v) for v in evr[:10]]
    res["pca_pc1_pct"] = float(evr[0] * 100)

    # ---------- 2. silhouette ----------
    sil = {}
    for k in range(2, 7):
        lab = KMeans(k, n_init=25, random_state=SEED).fit_predict(Z)
        sil[k] = float(silhouette_score(Z, lab))
    print("[2] silhouette（越接近 1 越好，<0.2 视为无明显结构）:",
          {k: round(v, 4) for k, v in sil.items()})
    res["silhouette"] = sil

    # ---------- 3. PAC vs 置换零分布 ----------
    obs = {}
    for k in [2, 3]:
        obs[k] = pac_score(M, k, 200, rng)
    null = {}
    for k in [2, 3]:
        vals = []
        for r in range(5):
            Mperm = rng.permuted(M, axis=0)  # 打散样本-事件配对（保留边缘分布）
            vals.append(pac_score(Mperm, k, 40, rng))
        null[k] = vals
    print("[3] PAC 观测值 vs 置换零分布:")
    for k in [2, 3]:
        mu, sd = np.mean(null[k]), np.std(null[k])
        z = (obs[k] - mu) / sd if sd > 0 else np.nan
        print(f"    k={k}: 观测 {obs[k]:.3f} | 零分布 {mu:.3f}±{sd:.3f} | z={z:+.2f}")
        res[f"pac_k{k}"] = obs[k]
        res[f"pac_null_k{k}_mean"] = float(mu)
        res[f"pac_null_k{k}_sd"] = float(sd)
        res[f"pac_k{k}_z"] = float(z) if sd > 0 else None

    # ---------- 4. RI 负荷分布形态 ----------
    b = samp["RI_burden"].values
    skew = float(st.skew(b))
    kurt = float(st.kurtosis(b))
    # 单峰性: 与同均值方差的正态分布比较 AIC（用 2/3 组分高斯混合的 BIC 差）
    from sklearn.mixture import GaussianMixture
    Xb = b.reshape(-1, 1)
    bic = {}
    for n in (1, 2, 3):
        gm = GaussianMixture(n, n_init=10, random_state=SEED).fit(Xb)
        bic[n] = float(gm.bic(Xb))
    print(f"\n[4] RI 负荷分布: 偏度={skew:.3f}, 峰度={kurt:.3f}")
    print(f"    GaussianMixture BIC: { {k: round(v,1) for k,v in bic.items()} }")
    print(f"    （BIC 最小者为最优组分数）")
    res["burden_skew"] = skew
    res["burden_kurtosis"] = kurt
    res["gmm_bic"] = bic
    res["burden_n_clusters_by_bic"] = int(min(bic, key=bic.get))
    # 置换检验：观测方差是否显著大于零假设
    print(f"    负荷范围: {b.min():.4f} – {b.max():.4f}（极差 {b.max()-b.min():.4f}）")
    res["burden_range"] = [float(b.min()), float(b.max())]

    # ---------- 5. 合并临床 / 生存 ----------
    clin = pd.read_csv(MATRIX, sep="\t", index_col=0, low_memory=False)
    clin.index = [barcode(i) for i in clin.index]
    clin = clin[~clin.index.duplicated()]
    surv = pd.read_csv(SURV, sep="\t", index_col=0)
    surv.index = [barcode(i) for i in surv.index]
    surv = surv[~surv.index.duplicated()]
    print(f"\n[5] clinicalMatrix: {clin.shape}  |  survival: {surv.shape}")

    keep = ["age_at_initial_pathologic_diagnosis", "gender", "histological_type",
            "neoplasm_histologic_grade", "GeneExp_Subtype", "G_CIMP_STATUS",
            "_primary_disease"]
    keep = [c for c in keep if c in clin.columns]
    tab = samp.join(clin[keep], how="left")
    tab = tab.join(surv[["OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time",
                         "PFI", "PFI.time"]], how="left")
    print(f"    合并后: {tab.shape}, 有 OS 数据者 {tab['OS'].notna().sum()}")
    n_both = int((tab["RI_group"].notna() & tab["OS"].notna()).sum())
    print(f"    RI 分组 + OS 均可用: {n_both} 例")
    res["n_with_os"] = n_both
    tab.to_csv(f"{OUT}/output/sample_RI_group_clinical.tsv", sep="\t")

    # ---------- 6. 与分子亚型的关联 ----------
    assoc = {}
    for col in ["GeneExp_Subtype", "G_CIMP_STATUS", "gender",
                "neoplasm_histologic_grade"]:
        if col not in tab.columns:
            continue
        sub = tab[[col, "RI_group"]].dropna()
        if sub[col].nunique() < 2:
            continue
        ct = pd.crosstab(sub["RI_group"], sub[col])
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            continue
        chi2, p, dof, _ = st.chi2_contingency(ct)
        assoc[col] = {"p": float(p), "chi2": float(chi2),
                      "table": ct.to_dict(), "n": int(len(sub))}
        print(f"\n[6] {col} × RI_group: p={p:.4f} (n={len(sub)})")
        print(ct.to_string())
    res["clinical_assoc"] = assoc

    # 年龄差异
    age = pd.to_numeric(tab["age_at_initial_pathologic_diagnosis"], errors="coerce")
    ok = age.notna() & tab["RI_group"].notna()
    if ok.sum() > 10:
        a = age[ok & (tab["RI_group"] == "RI-high")]
        c = age[ok & (tab["RI_group"] == "RI-low")]
        u, p = st.mannwhitneyu(a, c, alternative="two-sided")
        print(f"\n[6] 年龄: RI-high {a.mean():.1f} vs RI-low {c.mean():.1f}, p={p:.4f}")
        res["age_p"] = float(p)

    with open(f"{OUT}/output/step2a_structure.json", "w") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("\n完成 -> step2a_structure.json / sample_RI_group_clinical.tsv")


if __name__ == "__main__":
    main()
