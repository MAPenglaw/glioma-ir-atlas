#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""55 — 第二套 PSI 量化管道复现：TCGASpliceSeq (Ryan 2016) vs SplAdder (Kahles 2018)

【为什么要做】
B 线手稿 Limitations 第十三条自陈硬伤：全部 IR 结论建立在**单一 PSI 量化**
（Kahles 2018 pancan SplAdder）之上。本地已有第二种独立量化 ——
**TCGASpliceSeq**（Ryan et al. 2016, hg19, SpliceGraph 模型），数据在主线 A 项目，
零下载。

【两层 orthogonal 差异】
    SplAdder (Kahles 2018)  自建 splice graph，4 段坐标事件 ID，自 15 万+ IR 事件筛选
    SpliceSeq (Ryan 2016)   独立建图，as_id + exon 位点串，仅约 1.3k–2.3k 条 RI 事件
事件发现、内含子保留定义、定量读数全部不同 —— 唯一共享的是上游 TCGA RNA-seq。
所以本分析验证的是**定量阶段的管道稳健性**，不是"从 raw reads 重新定量"。

【矩阵口径严格对齐 50_residual_robustness.py / 16_km_core.py】
这是本脚本踩过的坑，务必保留：
  1. MIN_DETECT=0.90、MIN_STD=0.10，取高变异事件
  2. **逐事件 z-score** 后 nan_to_num(0)，转置成 样本 × 事件
     —— 漏掉 z-score 会让 silhouette 从 0.148 抬到 0.425，与手稿完全不可比
  3. silhouette 在**全高维空间**算，不在 PCA 投影空间算
  4. GaussianMixture 用 covariance_type="full"，不是 "diag"
  5. PAC 用修正版（分母=样本对实际共现次数），行列双 80% 重采样

【事件数不匹配的后果】两套管道高变异事件数差 3–7 倍，而 silhouette/PAC/BIC
对「样本数 : 特征数」比极度敏感。因此本脚本**自带下采样对照**：
把 SplAdder 的事件随机抽到 SpliceSeq 同等数量，重复多次，
判断差异究竟来自"管道"还是来自"特征数"。这是能否写进结论的前提。

用法  python 55_second_pipeline_spliceseq.py
输出  ri_atlas/output_secondpipeline/
"""
import os
import time

import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

ROOT_A = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc-as-glioma"
ROOT_B = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
ATLAS = os.path.join(ROOT_B, "ri_atlas")
OUT = os.path.join(ATLAS, "output_secondpipeline")
os.makedirs(OUT, exist_ok=True)

MIN_DETECT = 0.90
MIN_STD = 0.10
SEED = 20260917
SUB_SAMPLE = 0.80
SUB_FEATURE = 0.80
N_ITER_PAC = 200
N_DOWN_REP = 10

SPLICESEQ_SRC = {
    "GBM": ([f"{ROOT_A}/data/psi_batches/GBM/{i:03d}.txt" for i in range(19)]
            + [f"{ROOT_A}/data/psi/PSI_GBM_salvage.txt"]),
    "LGG": ([f"{ROOT_A}/data/psi_batches/LGG/{i:03d}.txt" for i in range(24)]
            + [f"{ROOT_A}/data/PSI_LGG_all.txt"]),
}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------- 矩阵
def to_X_aligned(vals):
    """严格复制 50_residual_robustness.load_hivar：
    逐事件 z-score → nan_to_num(0) → 转置为 样本 × 事件。"""
    vals = vals.apply(pd.to_numeric, errors="coerce")
    detect = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    M = vals.loc[(detect >= MIN_DETECT) & (sd >= MIN_STD)]
    Mv = M.values.astype(float)
    Xz = (Mv - np.nanmean(Mv, axis=1, keepdims=True)) / (np.nanstd(Mv, axis=1, keepdims=True) + 1e-12)
    X = np.nan_to_num(Xz, nan=0.0).T            # 样本 × 事件
    return X, np.array([str(c) for c in M.columns])


def load_spliceseq(cancer):
    frames = []
    for p in SPLICESEQ_SRC[cancer]:
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t", low_memory=False)
            frames.append(d[d["splice_type"] == "RI"])
    m = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["as_id"])
    vals = m.drop(columns=["symbol", "as_id", "splice_type", "exons", "from_exon",
                           "to_exon", "novel_splice", "pct_with_values", "psi_range",
                           "std_psi"], errors="ignore")
    vals.columns = [str(c).replace("_", "-")[:12] for c in vals.columns]
    vals = vals.T.groupby(level=0).mean().T
    return vals


def load_spladder(cancer):
    psi = pd.read_csv(f"{ATLAS}/data/KM_IR_{cancer}.psi.gz", sep="\t",
                      compression="gzip", low_memory=False)
    vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    vals = vals.loc[:, ~vals.columns.duplicated()]
    vals.columns = [str(c).split(".")[0].strip().strip('"')[:12] for c in vals.columns]
    return vals.T.groupby(level=0).mean().T


def load_lineb_burden(cancer):
    b = pd.read_csv(f"{ATLAS}/output_km/{cancer}/burden.tsv", sep="\t", index_col=0)
    b.index = [str(i)[:12] for i in b.index]
    return b["RI_burden"].groupby(level=0).mean()


# ---------------------------------------------------------------- 诊断
def pac_k2(X, n_iter=N_ITER_PAC, seed=SEED):
    """修正版 PAC（Monti 2003 分母），行列双 80% 重采样。"""
    rng = np.random.default_rng(seed)
    n_s, p = X.shape
    if n_s < 10 or p < 2:
        return float("nan")
    n_sub = max(3, int(round(SUB_SAMPLE * n_s)))
    p_sub = max(2, int(round(SUB_FEATURE * p)))
    conn = np.zeros((n_s, n_s))
    valid = np.zeros((n_s, n_s))
    for it in range(n_iter):
        rs = rng.choice(n_s, n_sub, replace=False)
        cf = rng.choice(p, p_sub, replace=False)
        sub = X[np.ix_(rs, cf)]
        lab = KMeans(n_clusters=2, n_init=10, random_state=seed + it).fit_predict(sub)
        pos = {s: i for i, s in enumerate(rs)}
        for a in range(len(rs)):
            for b in range(a + 1, len(rs)):
                sa, sb = rs[a], rs[b]
                same = 1 if lab[pos[sa]] == lab[pos[sb]] else 0
                conn[sa, sb] += same
                conn[sb, sa] += same
                valid[sa, sb] += 1
                valid[sb, sa] += 1
    C = np.where(valid > 0, conn / np.maximum(valid, 1), 0.0)
    np.fill_diagonal(C, 1.0)
    v = C[np.triu_indices(n_s, k=1)]
    return float(np.mean((v > 0.1) & (v < 0.9)))


def diagnostics(X, burden, tag, pac=True):
    """严格对齐 50_residual_robustness.diagnostics()。"""
    n_s = X.shape[0]
    ncomp = min(10, X.shape[0] - 1, X.shape[1])
    pc = PCA(n_components=ncomp, random_state=SEED).fit_transform(X)
    pc1 = pc[:, 0]
    if st.spearmanr(pc1, burden)[0] < 0:
        pc1 = -pc1
    evr = float(PCA(n_components=ncomp, random_state=SEED).fit(X).explained_variance_ratio_[0])
    lab = KMeans(n_clusters=2, n_init=50, random_state=SEED).fit_predict(X)
    sil = float(silhouette_score(X, lab, sample_size=min(n_s, 1000)))
    bic = {}
    for k in (1, 2):
        bic[k] = float(GaussianMixture(n_components=k, covariance_type="full",
                                       random_state=SEED).fit(X).bic(X))
    return {"pipe": tag, "n_sample": int(n_s), "n_event": int(X.shape[1]),
            "PC1_var": round(evr, 4),
            "rho_PC1_burden": round(float(st.spearmanr(pc1, burden)[0]), 4),
            "silhouette_k2": round(sil, 4),
            "BIC_k1": round(bic[1], 1), "BIC_k2": round(bic[2], 1),
            "BIC_prefers_k1": bool(bic[1] < bic[2]),
            "PAC_k2": round(pac_k2(X), 4) if pac else np.nan}


# ---------------------------------------------------------------- 主流程
def main():
    summary, pairs = [], []
    for cancer in ["GBM", "LGG"]:
        log(f"\n{'=' * 70}\n{cancer}\n{'=' * 70}")

        # ---- 第一套：SplAdder ----
        va = load_spladder(cancer)
        ba_line = load_lineb_burden(cancer)
        common = sorted(set(va.columns) & set(ba_line.index))
        va = va[common]
        Xa, _ = to_X_aligned(va)
        burden_a = va.loc[(va.notna().mean(axis=1) >= MIN_DETECT)].mean(axis=0, skipna=True).values
        log(f"  SplAdder : 事件池 {va.shape[0]} → 高变异 {Xa.shape[1]}，样本 {Xa.shape[0]}")
        da = diagnostics(Xa, burden_a, "SplAdder")
        log(f"    PC1={da['PC1_var']:.3f} sil={da['silhouette_k2']:.3f} "
            f"PAC={da['PAC_k2']:.3f} BIC偏好={'k1' if da['BIC_prefers_k1'] else 'k2'}")

        # ---- 第二套：SpliceSeq ----
        vs = load_spliceseq(cancer)
        vs = vs[[c for c in common if c in vs.columns]]
        common2 = sorted(set(vs.columns) & set(common))
        va2 = va[common2]
        vs = vs[common2]
        Xa2, _ = to_X_aligned(va2)
        burden_a2 = va2.loc[(va2.notna().mean(axis=1) >= MIN_DETECT)].mean(axis=0, skipna=True).values
        Xs, _ = to_X_aligned(vs)
        burden_s = vs.loc[(vs.notna().mean(axis=1) >= MIN_DETECT)].mean(axis=0, skipna=True).values
        log(f"  SpliceSeq: 事件池 {vs.shape[0]} → 高变异 {Xs.shape[1]}，样本 {Xs.shape[0]}")
        ds = diagnostics(Xs, burden_s, "SpliceSeq")
        log(f"    PC1={ds['PC1_var']:.3f} sil={ds['silhouette_k2']:.3f} "
            f"PAC={ds['PAC_k2']:.3f} BIC偏好={'k1' if ds['BIC_prefers_k1'] else 'k2'}")

        # ---- T1 两套 burden 一致性（同一批共同患者）----
        b_ss = pd.Series(burden_s, index=common2)
        b_line = ba_line.reindex(common2)
        rho, pv = st.spearmanr(b_ss.values, b_line.values)
        r_p, pv_p = st.pearsonr(b_ss.values, b_line.values)
        log(f"  [T1] 两套 burden：Spearman={rho:.3f} (p={pv:.2e})  Pearson={r_p:.3f} (n={len(common2)})")

        pairs.append(pd.DataFrame({"cancer": cancer, "sample": common2,
                                   "burden_SplAdder": b_line.values,
                                   "burden_SpliceSeq": b_ss.values}))

        # ---- 下采样对照：SplAdder 事件抽到 SpliceSeq 同等数量 ----
        tgt = min(Xs.shape[1], Xa2.shape[1])
        downs = []
        for rep in range(N_DOWN_REP):
            cols = np.random.default_rng(SEED + rep).choice(Xa2.shape[1], size=tgt, replace=False)
            dd = diagnostics(Xa2[:, cols], burden_a2, f"SplAdder_ds{rep}", pac=False)
            downs.append(dd)
            if (rep + 1) % 5 == 0:
                log(f"    下采样 {rep + 1}/{N_DOWN_REP}: sil={dd['silhouette_k2']:.3f} "
                    f"BIC={'k1' if dd['BIC_prefers_k1'] else 'k2'}")
        sil_arr = np.array([d["silhouette_k2"] for d in downs])
        k2 = sum(1 for d in downs if not d["BIC_prefers_k1"])
        log(f"    下采样汇总 sil={sil_arr.mean():.3f}±{sil_arr.std():.3f}，BIC偏好k2比例={k2}/{len(downs)}")

        summary.append({
            "cancer": cancer, "n_common": len(common2),
            "T1_burden_spearman": round(float(rho), 3),
            "T1_burden_p": f"{pv:.1e}",
            "SA_event": da["n_event"], "SS_event": ds["n_event"],
            "SA_PC1": da["PC1_var"], "SS_PC1": ds["PC1_var"],
            "SA_rho_PC1_burden": da["rho_PC1_burden"], "SS_rho_PC1_burden": ds["rho_PC1_burden"],
            "SA_sil": da["silhouette_k2"], "SS_sil": ds["silhouette_k2"],
            "SA_ds_sil_mean": round(float(sil_arr.mean()), 3),
            "SA_ds_sil_sd": round(float(sil_arr.std()), 3),
            "SA_PAC": da["PAC_k2"], "SS_PAC": ds["PAC_k2"],
            "SA_BIC_k1": da["BIC_prefers_k1"], "SS_BIC_k1": ds["BIC_prefers_k1"],
            "SA_ds_BIC_k2_frac": round(k2 / len(downs), 2),
        })

    R = pd.DataFrame(summary)
    R.to_csv(f"{OUT}/cross_pipeline_results.csv", index=False)
    pd.concat(pairs, ignore_index=True).to_csv(f"{OUT}/burden_crosscheck.csv", index=False)
    log("\n" + R.to_string(index=False))
    log(f"\n输出 {OUT}/cross_pipeline_results.csv")


if __name__ == "__main__":
    main()
