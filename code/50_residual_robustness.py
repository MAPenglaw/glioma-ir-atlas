#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""50 — 残差化稳健性检验（投前必做，评审 P0-1）。

问题：手稿披露了 burden~TotalExpr (GBM +0.459 / LGG +0.298)，但**未披露更强的**
PC1~TotalExpr (GBM +0.551 / LGG +0.348)。审稿人重算即可发现。
本脚本回答唯一可能动摇地基的问题：

    **把总表达量从 PSI 矩阵里彻底扣掉之后，"连续轴"还在不在？**

做法：对高变异 IR 事件矩阵 X (样本 × 事件，逐事件 z-score)，逐事件对混杂做 OLS
取残差，得到 X_res，然后**重跑全部轴诊断**，与原始口径对照。

三个口径：
  R0  原始（对照，复现既有数字）
  R1  扣掉 TotalExpr（全基因 log2(TPM+0.001) 均值，技术混杂主代理）
  R2  扣掉 TotalExpr + Tumor_content + Neuronal_normal_brain（最严格）

诊断指标（每口径）：
  PC1 方差解释率、PC1~burden、PC1~TotalExpr、PC1~Sm_core、
  silhouette(k=2)、PAC(k=2)、GMM BIC(k=1 vs k=2)

输出：
  output_km/residual_robustness.tsv    长表 metric × cancer × condition
  output_km/residual_robustness.json   同上 + 判定
  output_km/figS3_residual_scores.tsv  补充图数据（PC1 / burden / TotalExpr）
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
DATA = os.path.join(ROOT, "ri_atlas", "data")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
CONF = os.path.join(ROOT, "ri_atlas", "output_pancan", "tme_confounders.tsv")

MIN_DETECT = 0.90
MIN_STD = 0.10
N_ITER = 100          # consensus 重采样次数（与 22 脚本一致）
SUB_SAMPLE = 0.8
SUB_FEATURE = 0.8
SEED = 20260927


# ---------------------------------------------------------------- 数据
def load_hivar(cancer):
    """返回 X (样本 × 高变异事件，逐事件 z-score) 与样本 barcode。"""
    psi = pd.read_csv(os.path.join(DATA, f"KM_IR_{cancer}.psi.gz"), sep="\t",
                      compression="gzip", low_memory=False)
    vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    vals = vals.loc[:, ~vals.columns.duplicated()]
    detect = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    mask = (detect >= MIN_DETECT) & (sd >= MIN_STD)
    M = vals.loc[mask]
    Mv = M.values.astype(float)
    mu = np.nanmean(Mv, axis=1, keepdims=True)
    sdd = np.nanstd(Mv, axis=1, keepdims=True)
    Xz = (Mv - mu) / (sdd + 1e-12)
    X = np.nan_to_num(Xz, nan=0.0).T                      # 样本 × 事件
    cols = [str(c).strip().strip('"')[:12] for c in vals.columns]
    return X, np.array(cols), int(mask.sum())


def ranknorm(v):
    """秩 → 标准正态分位（单调不变，抗离群）。"""
    v = np.asarray(v, dtype=float)
    r = st.rankdata(v)
    return st.norm.ppf(r / (len(r) + 1.0))


def residualize(X, C):
    """逐事件对设计矩阵 C (样本 × 协变量，含截距) 取 OLS 残差，再重新 z-score。"""
    C = np.column_stack([np.ones(X.shape[0]), C])
    # 最小二乘：beta = (C'C)^-1 C' X
    beta, *_ = np.linalg.lstsq(C, X, rcond=None)
    Xr = X - C @ beta
    mu = Xr.mean(axis=0, keepdims=True)
    sd = Xr.std(axis=0, keepdims=True)
    return (Xr - mu) / (sd + 1e-12)


# ---------------------------------------------------------------- 诊断
def pac_k2(X):
    """consensus clustering PAC（k=2），与 22 脚本同一算法。"""
    rng = np.random.default_rng(SEED)
    n_s, n_f = X.shape
    conn = np.zeros((n_s, n_s))
    valid = np.zeros((n_s, n_s))
    for it in range(N_ITER):
        rs = rng.choice(n_s, int(n_s * SUB_SAMPLE), replace=False)
        fs = rng.choice(n_f, int(n_f * SUB_FEATURE), replace=False)
        sub = X[np.ix_(rs, fs)]
        km = KMeans(n_clusters=2, n_init=10, random_state=SEED + it)
        lab = km.fit_predict(sub)
        idx = np.where(np.isin(np.arange(n_s), rs))[0]
        pos = {s: i for i, s in enumerate(rs)}
        for a in range(len(rs)):
            for b in range(a + 1, len(rs)):
                sa, sb = rs[a], rs[b]
                conn[sa, sb] += (lab[a] == lab[b])
                conn[sb, sa] += (lab[a] == lab[b])
                valid[sa, sb] += 1
                valid[sb, sa] += 1
    C = np.where(valid > 0, conn / np.maximum(valid, 1), 0.0)
    np.fill_diagonal(C, 1.0)
    iu = np.triu_indices(n_s, k=1)
    v = C[iu]
    return float(np.mean((v > 0.1) & (v < 0.9)))


def diagnostics(X, burden, tot, smc):
    n_s = X.shape[0]
    ncomp = min(10, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=ncomp)
    pc = pca.fit_transform(X)
    pc1 = pc[:, 0]
    evr = float(pca.explained_variance_ratio_[0])
    # PC1 方向取与 burden 正相关（保证符号可比）
    if st.spearmanr(pc1, burden)[0] < 0:
        pc1 = -pc1

    km2 = KMeans(n_clusters=2, n_init=50, random_state=SEED)
    lab = km2.fit_predict(X)
    sil = float(silhouette_score(X, lab, sample_size=min(n_s, 1000)))
    bic = {}
    for k in (1, 2):
        g = GaussianMixture(n_components=k, covariance_type="full", random_state=SEED)
        g.fit(X)
        bic[k] = float(g.bic(X))

    def sp(a, b):
        if a is None:
            return (np.nan, np.nan)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 10:
            return (np.nan, np.nan)
        r, p = st.spearmanr(a[m], b[m])
        return (float(r), float(p))

    r_b, p_b = sp(pc1, burden)
    r_t, p_t = sp(pc1, tot)
    r_s, p_s = sp(pc1, smc)
    r_bt, p_bt = sp(burden, tot)
    return {
        "n_samples": n_s,
        "PC1_var": evr,
        "PC2_var": float(pca.explained_variance_ratio_[1]),
        "silhouette_k2": sil,
        "PAC_k2": pac_k2(X),
        "BIC_k1": bic[1],
        "BIC_k2": bic[2],
        "BIC_prefers_k1": bool(bic[1] < bic[2]),
        "rho_PC1_burden": r_b, "p_PC1_burden": p_b,
        "rho_PC1_TotalExpr": r_t, "p_PC1_TotalExpr": p_t,
        "rho_PC1_Sm_core": r_s, "p_PC1_Sm_core": p_s,
        "rho_burden_TotalExpr": r_bt, "p_burden_TotalExpr": p_bt,
        "_pc1": pc1,
    }


# ---------------------------------------------------------------- 主流程
def main():
    os.makedirs(OUT, exist_ok=True)
    conf = pd.read_csv(CONF, sep="\t")
    conf["key"] = conf["sample"].astype(str).str[:12]
    conf = conf.drop_duplicates("key").set_index("key")

    rows, plot_rows = [], []
    for cancer in ("GBM", "LGG"):
        X, barcodes, n_ev = load_hivar(cancer)
        print(f"[{cancer}] 高变异事件 {n_ev}，样本 {X.shape[0]}", flush=True)

        cf = conf.reindex(barcodes)
        ok = cf["RI_burden"].notna().values
        X = X[ok]
        barcodes = barcodes[ok]
        cf = cf.loc[cf["RI_burden"].notna()]
        cf = cf.reindex(barcodes)
        print(f"  匹配混杂变量 {X.shape[0]} 样本", flush=True)

        burden = cf["RI_burden"].values.astype(float)
        tot = cf["TotalExpr"].values.astype(float)
        smc = cf["Sm_core"].values.astype(float)
        purity = cf["Tumor_content"].values.astype(float)
        neuron = cf["Neuronal_normal_brain"].values.astype(float)

        conds = {
            "R0_raw": X,
            "R1_resid_TotalExpr": residualize(X, ranknorm(tot).reshape(-1, 1)),
            "R2_resid_full": residualize(
                X, np.column_stack([ranknorm(tot), ranknorm(purity), ranknorm(neuron)])),
        }

        res = {}
        for name, Xc in conds.items():
            print(f"  --- {name} ---", flush=True)
            d = diagnostics(Xc, burden, tot, smc)
            res[name] = d
            for kk in ("PC1_var", "PC2_var", "silhouette_k2", "PAC_k2",
                       "BIC_k1", "BIC_k2", "BIC_prefers_k1",
                       "rho_PC1_burden", "p_PC1_burden",
                       "rho_PC1_TotalExpr", "p_PC1_TotalExpr",
                       "rho_PC1_Sm_core", "p_PC1_Sm_core",
                       "rho_burden_TotalExpr", "p_burden_TotalExpr"):
                rows.append({"cancer": cancer, "condition": name, "metric": kk,
                             "value": d[kk], "n_samples": d["n_samples"],
                             "n_events": n_ev})
            plot_rows.append(pd.DataFrame({
                "sample": barcodes, "cancer": cancer, "condition": name,
                "PC1": d["_pc1"], "RI_burden": burden,
                "TotalExpr": tot, "Sm_core": smc}))

        # 关键判定
        r0, r1, r2 = res["R0_raw"], res["R1_resid_TotalExpr"], res["R2_resid_full"]
        verdict = {
            "cancer": cancer,
            "PC1_var_R0": r0["PC1_var"], "PC1_var_R1": r1["PC1_var"], "PC1_var_R2": r2["PC1_var"],
            "PC1_var_retained_frac_R2": r2["PC1_var"] / r0["PC1_var"],
            "PAC_R0": r0["PAC_k2"], "PAC_R1": r1["PAC_k2"], "PAC_R2": r2["PAC_k2"],
            "sil_R0": r0["silhouette_k2"], "sil_R1": r1["silhouette_k2"], "sil_R2": r2["silhouette_k2"],
            "BIC_k1_better_R2": r2["BIC_prefers_k1"],
            "rho_PC1_burden_R0": r0["rho_PC1_burden"], "rho_PC1_burden_R2": r2["rho_PC1_burden"],
            "rho_PC1_TotalExpr_R0": r0["rho_PC1_TotalExpr"],
            "rho_PC1_TotalExpr_R2": r2["rho_PC1_TotalExpr"],
            "rho_PC1_Sm_core_R0": r0["rho_PC1_Sm_core"],
            "rho_PC1_Sm_core_R2": r2["rho_PC1_Sm_core"],
            "p_PC1_Sm_core_R2": r2["p_PC1_Sm_core"],
        }
        # 判定：轴仍在 = 扣掉总表达后 PC1 方差解释率仍保留 ≥50% 且 PC1~burden 仍 ≥0.8
        verdict["AXIS_SURVIVES"] = bool(
            verdict["PC1_var_retained_frac_R2"] >= 0.5 and r2["rho_PC1_burden"] >= 0.80)
        # 判定：仍无离散亚型 = 残差化后 PAC 仍高（>0.4）且 BIC 仍偏好 k=1
        verdict["STILL_NO_DISCRETE_SUBTYPE"] = bool(
            r2["PAC_k2"] > 0.4 and r2["BIC_prefers_k1"])
        print(f"  >>> 轴存活={verdict['AXIS_SURVIVES']}  "
              f"PC1方差保留={verdict['PC1_var_retained_frac_R2']:.2%}  "
              f"R2 PC1~burden={r2['rho_PC1_burden']:+.3f}  "
              f"R2 PC1~TotalExpr={r2['rho_PC1_TotalExpr']:+.3f}  "
              f"R2 PC1~Sm_core={r2['rho_PC1_Sm_core']:+.3f} (p={r2['p_PC1_Sm_core']:.2e})",
              flush=True)
        res["_verdict"] = verdict

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "residual_robustness.tsv"),
                              sep="\t", index=False)
    pd.concat(plot_rows).to_csv(os.path.join(OUT, "figS3_residual_scores.tsv"),
                                sep="\t", index=False)
    # _pc1 是 numpy 数组，仅供绘图，不入 json
    json_safe = {c: {k: (v if k != "_pc1" else None)
                     for k, v in d.items()} for c, d in res.items()}
    with open(os.path.join(OUT, "residual_robustness.json"), "w") as f:
        json.dump(json_safe, f, indent=1, default=float)
    print("DONE -> output_km/residual_robustness.{tsv,json} + figS3_residual_scores.tsv")


if __name__ == "__main__":
    main()
