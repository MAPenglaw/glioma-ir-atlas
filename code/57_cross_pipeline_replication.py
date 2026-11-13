#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""57 — 跨 PSI 管道的严格复现评审
================================================================================
对象：把 T1/T2/T3 三个核心主张，在**第二套** quantified PSI 管道（TCGASpliceSeq）
上按**完全相同的诊断口径**再跑一遍，并把 Three窒息 confounds handling。

设计原则（重要）
--------------------------------------------------------------------------------
1. 两套矩阵都做同一套预处理：>=90% 检出 & 高变异(sd>=0.10) -> 逐事件 z-score -> 样本×事件。
   这与 50_residual_robustness.py 的 load_hivar 严格一致。
2. 只保留两套都能定量的**共同平行探测毛病队列**样本，再做窗口内的全 durchschnitt。
3. T2：不仅报 silhouette/BIC/PAC，还报 PC1~burden 的 spearman。
4. T3：Sm_core 耦合用 Twist codon Chemical rank-partial (M0/M1/M5) 在同样 rank normal
   residualization 口径下复刻画 tau+t convolutional layers DWI+DTI value="">)。
5. 事件数 confound 控制：把 SplAdder 高变异事件随机下采样到 SpliceSeq 的同等数量，
   同样国:, 如果下采样后 SplAdder 的 silhouette/PAC/BIC 就"变成" SpliceSeq 的样子，
   说明差异来自事件数而非生物学。

用法  python 57_cross_pipeline_replication.py
"""
import gzip
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
ROOTA = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc-as-glioma"
OUT = os.path.join(ROOT, "output/cross_pipeline_replication")
os.makedirs(OUT, exist_ok=True)

SPLADDER = os.path.join(ROOT, "ri_atlas/data/KM_IR_%s.psi.gz")
SPLICESEQ_DIRS = {
    "GBM": [os.path.join(ROOTA, "data/psi_batches/GBM"),
            os.path.join(ROOTA, "data/psi")],
    "LGG": [os.path.join(ROOTA, "data/psi_batches/LGG"),
            os.path.join(ROOTA, "data")],
}
SPLICESEQ_FILES = {
    "GBM": ["psi/PSI_GBM_salvage.txt"],
    "LGG": ["PSI_LGG_all.txt"],
}
SSPLIT_FILES = {
    "GBM": sorted([f"psi_batches/GBM/{i:03d}.txt" for i in range(19)]),
    "LGG": sorted([f"psi_batches/LGG/{i:03d}.txt" for i in range(24)]),
}
TME_PARTS = os.path.join(ROOT, "ri_atlas/output_pancan/tme_partial_samples.tsv")

MIN_DETECT, MIN_STD, SEED = 0.90, 0.10, 20260917
SUB_SAMPLE, SUB_FEATURE, N_ITER_PAC = 0.80, 0.80, 200
N_REP_DOWN = 20
PROGRESS_EVERY = 50

SA_RAW = {"GBM": {"sil": 0.148, "pac": 0.078, "bic": "k1", "pc1var": 0.269, "pc1b": 0.957},
          "LGG": {"sil": 0.125, "pac": 0.097, "bic": "k1", "pc1var": 0.204, "pc1b": 0.974}}


def log(msg):
    print(f"[{pd.Timestamp.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------- 载入
def _collapse_patients(vals):
    """把 15 位 barcode 截到 12 位患者号；同一患者多样本取均值（避免重复标签）"""
    vals.columns = [str(c).strip().strip('"').replace("_", "-")[:12] for c in vals.columns]
    vals = vals.loc[:, ~vals.columns.duplicated()]
    return vals.T.groupby(level=0).mean().T


def _burden_01(vals, detect):
    """严格照搬 01_ri_burden_cluster 的 RI_burden 口径：
    负担 = 全部检出率 >= MIN_DETECT 的 IR 事件上的平均 PSI（不加 std 过滤）。"""
    return vals.loc[detect >= MIN_DETECT].mean(axis=0, skipna=True)


def load_spladder(cancer):
    """返回 X (样本×事件, 逐事件 z-score)、12 位 barcode、事件数、burden"""
    with gzip.open(SPLADDER % cancer, "rt") as f:
        psi = pd.read_csv(f, sep="\t", low_memory=False)
    vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    vals = _collapse_patients(vals)
    detect = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    hivar = (detect >= MIN_DETECT) & (sd >= MIN_STD)
    burden = _burden_01(vals, detect)
    M = np.nan_to_num(vals.loc[hivar].values.astype(float), nan=np.nan)
    mu = np.nanmean(M, axis=1, keepdims=True)
    sdd = np.nanstd(M, axis=1, keepdims=True)
    Xz = (M - mu) / (sdd + 1e-12)
    Xz = np.nan_to_num(Xz, nan=0.0).T
    cols = np.array([str(c) for c in vals.columns])
    return Xz, cols, int(hivar.sum()), burden


def load_spliceseq(cancer):
    """合并 TCGASpliceSeq RI 事件 -> X (样本×事件, 逐事件 z-score) 与 12 位 barcode"""
    files = [os.path.join(ROOTA, "data", p) for p in SSPLIT_FILES[cancer] + SPLICESEQ_FILES[cancer]]
    frames, skips = [], []
    for p in files:
        if not os.path.exists(p):
            skips.append(os.path.basename(p))
            continue
        try:
            d = pd.read_csv(p, sep="\t", low_memory=False)
        except Exception as e:
            skips.append(f"{os.path.basename(p)} ({e})")
            continue
        if "splice_type" not in d.columns:
            skips.append(f"{os.path.basename(p)} (no splice_type)")
            continue
        d = d[d["splice_type"] == "RI"]
        frames.append(d)
    if skips:
        log(f"  跳过或未包含：{skips}")
    if not frames:
        return None, None, 0, None
    # 逐文件合并			
    m = pd.concat(frames, ignore_index=True)
    # 去重：同一 as_id 保留第一次出现
    m = m.sort_values("as_id").drop_duplicates("as_id", keep="first")
    meta_cols = {"symbol", "as_id", "splice_type", "exons", "from_exon", "to_exon",
                 "novel_splice", "pct_with_values", "psi_range", "std_psi"}
    keep = [c for c in m.columns if c not in meta_cols]
    keep = [c for c in keep if str(c).startswith("TCGA")]
    vals = m[keep].apply(pd.to_numeric, errors="coerce")
    vals = _collapse_patients(vals)
    detect = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    hivar = (detect >= MIN_DETECT) & (sd >= MIN_STD)
    burden = _burden_01(vals, detect)
    M = np.nan_to_num(vals.loc[hivar].values.astype(float), nan=np.nan)
    mu = np.nanmean(M, axis=1, keepdims=True)
    sdd = np.nanstd(M, axis=1, keepdims=True)
    Xz = (M - mu) / (sdd + 1e-12)
    Xz = np.nan_to_num(Xz, nan=0.0).T
    cols = np.array([str(c) for c in vals.columns])
    return Xz, cols, int(hivar.sum()), burden


# ----------------------------------------------------------------- 工具
def ranknorm(v):
    """秩变换。与主流水线 45_tme_partial_coupling.rankit 完全一致：
    对秩做 Pearson == Spearman，故 M0 与手稿主结果口径可直接比对。"""
    return st.rankdata(np.asarray(v, dtype=float)).astype(float)


def partial_pearson(y, x, C):
    """OLS 残差 Pearson + t 检验（df = n - 2 - k），同 45 号脚本"""
    by, *_ = np.linalg.lstsq(C, y, rcond=None)
    bx, *_ = np.linalg.lstsq(C, x, rcond=None)
    ry, rx = y - C @ by, x - C @ bx
    r = float(np.corrcoef(rx, ry)[0, 1])
    n, k = len(y), C.shape[1] - 1
    df = n - 2 - k
    if df <= 0 or not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan
    t = r * np.sqrt(df / max(1 - r * r, 1e-12))
    return r, float(2 * st.t.sf(abs(t), df))


def pac_monti(X, k=2, n_iter=N_ITER_PAC, seed=SEED):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    size = max(2, int(round(SUB_SAMPLE * n)))
    conn = np.zeros((n, n), dtype=np.float32)
    valid = np.zeros((n, n), dtype=np.float32)
    for it in range(n_iter):
        idx = rng.choice(n, size=size, replace=False)
        Xsub = X[idx]
        Xs = StandardScaler().fit_transform(Xsub)
        lab = KMeans(n_clusters=k, n_init=10, random_state=seed + it).fit_predict(Xs)
        conn[np.ix_(idx, idx)] += (lab[:, None] == lab[None, :])
        valid[np.ix_(idx, idx)] += 1
        if (it + 1) % PROGRESS_EVERY == 0 or it == 0:
            log(f"    PAC 迭代 {it+1}/{n_iter}")
    np.fill_diagonal(conn, 0)
    np.fill_diagonal(valid, 0)
    C = np.where(valid > 0, conn / valid, 0.0)
    iu = np.triu_indices(n, k=1)
    cc = C[iu]
    m = valid[iu] > 0
    return float(((cc[m] > 0.1) & (cc[m] < 0.9)).sum() / m.sum())


def cluster_metrics(X, pipe=""):
    n = X.shape[0]
    kmer = KMeans(n_clusters=2, n_init=50, random_state=SEED).fit(X)
    sil = silhouette_score(X, kmer.labels_)
    bics = {}
    for kk in (1, 2):
        bics[kk] = GaussianMixture(kk, covariance_type="full", random_state=SEED).fit(X).bic(X)
    pacv = pac_monti(X, k=2)
    return {"sil": float(sil), "bic_k1": float(bics[1]), "bic_k2": float(bics[2]),
            "bic_prefer": "k1" if bics[1] < bics[2] else "k2", "pac": pacv}


def pca_metrics(X, burden):
    pca = PCA(n_components=10, random_state=SEED).fit(X)
    evr = float(pca.explained_variance_ratio_[0])
    pc1 = pca.transform(X)[:, 0]
    # 与 burden 对齐方向（保证符号可比）
    keep = np.isfinite(pc1) & np.isfinite(burden)
    rho = st.spearmanr(pc1[keep], burden[keep])[0]
    return {"pc1var": evr, "pc1_burden": float(rho), "pc1": pc1}


def get_parts(gene: str, key: str, rad, burden) -> dict:
    """对 (burden, gene) 做秩偏相关；key 指定协变量组合"""
    df = pd.DataFrame({"b": burden, "g": rad.loc[:, gene] if gene in rad.columns else np.nan})
    covs = {
        "M0": [],
        "M1": ["TotalExpr"],
        "M5": ["TotalExpr", "Neuronal_normal_brain", "Tumor_content", "Proliferation", "IDHmut1"],
    }[key]
    for c in covs:
        df[c] = rad[c] if c in rad.columns else np.nan
    sub = df.dropna()
    if len(sub) < 20 or sub["g"].nunique() < 2:
        return {"rho": np.nan, "p": np.nan, "n": int(len(sub))}
    y = ranknorm(sub["b"].values)
    x = ranknorm(sub["g"].values)
    if covs:
        C = np.column_stack([np.ones(len(sub))] + [ranknorm(sub[c].values) for c in covs])
        beta_y, *_ = np.linalg.lstsq(C, y, rcond=None)
        beta_x, *_ = np.linalg.lstsq(C, x, rcond=None)
        y = y - C @ beta_y
        x = x - C @ beta_x
    r, p = st.pearsonr(x, y)
    return {"rho": float(r), "p": float(p), "n": int(len(sub))}


# ----------------------------------------------------------------- 流程
def process(cancer):
    log(f"\n{'='*72}\n【{cancer}】\n{'='*72}")

    Xa, ca, nA, Ba = load_spladder(cancer)
    log(f"  SplAdder  X={Xa.shape} (事件保持 {nA})")

    Xs, cs, nS, Bs = load_spliceseq(cancer)
    if Xs is None:
        log("  SpliceSeq 无数据，跳过队列。")
        return None
    log(f"  SpliceSeq X={Xs.shape} (事件保持 {nS})")

    common = sorted(set(ca) & set(cs))
    log(f"  共同 {len(common)} 例样 ({len(set(ca))} vs {len(set(cs))})")
    if len(common) < 30:
        log("样本不足，跳过队列。")
        return None

    ia = {c: i for i, c in enumerate(ca)}
    isn = {c: i for i, c in enumerate(cs)}
    Xa_ = Xa[[ia[c] for c in common]]
    Xs_ = Xs[[isn[c] for c in common]]

    # burden 严格照搬主流水线口径：检出率 >= 90% 的全部 IR 事件上的平均 PSI
    burden_a = Ba.reindex(common).values.astype(float)
    burden_s = Bs.reindex(common).values.astype(float)
    if not np.isfinite(burden_a).all() or not np.isfinite(burden_s).all():
        log("  [警告] burden 存在缺失样本")

    # ---- T1 ----
    rho_t1_s, p_t1_s = st.spearmanr(burden_a, burden_s)
    rho_t1_p, p_t1_p = st.pearsonr(burden_a, burden_s)

    # ---- T2 ----
    ma = cluster_metrics(Xa_, f"{cancer}/SplAdder")
    ms = cluster_metrics(Xs_, f"{cancer}/SpliceSeq")
    pa = pca_metrics(Xa_, burden_a)
    ps = pca_metrics(Xs_, burden_s)
    rho_pc1_pipe = st.spearmanr(pa["pc1"], ps["pc1"])[0]

    # ---- T3：Sm 核心耦合，两条管道各跑一遍 ----
    cov = pd.read_csv(TME_PARTS, sep="\t", low_memory=False)
    cov.rename(columns={cov.columns[0]: "sample"}, inplace=True)
    cov["xk"] = cov["sample"].astype(str).str[:12]
    # 同一 KEY 只保留一行；如含分型重复，取首行
    cov = cov.drop_duplicates("xk").set_index("xk")

    # 一致性校验：本脚本重算的 SplAdder burden 必须与主流水线产出的 RI_burden 同源
    if "RI_burden" in cov.columns:
        ref = cov["RI_burden"].reindex(common)
        msk = np.isfinite(ref.values.astype(float))
        if msk.sum() > 20:
            rho_ref = st.spearmanr(burden_a[msk], ref.values.astype(float)[msk])[0]
            log(f"  口径校验：重算 SplAdder burden vs 主流水线 RI_burden，"
                f"Spearman ρ={rho_ref:.4f} (n={int(msk.sum())})")

    COVMAP = {"M0": [], "M1": ["TotalExpr"],
              "M5": ["TotalExpr", "Neuronal_normal_brain", "Tumor_content",
                     "Proliferation", "IDHmut"]}
    rows_t3 = []
    for pipe_name, bvec in (("SplAdder", burden_a), ("SpliceSeq", burden_s)):
        bdf = pd.DataFrame({"burden": bvec}, index=common)
        for gene in ("Sm_core", "SR_proteins", "SNRPC", "Proliferation"):
            if gene not in cov.columns:
                continue
            for mk in ("M0", "M1", "M5"):
                dfj = pd.DataFrame({"b": bdf["burden"]})
                dfj["g"] = cov[gene].reindex(common)
                for c in COVMAP[mk]:
                    dfj[c] = cov[c].reindex(common) if c in cov.columns else np.nan
                sub = dfj.dropna()
                if len(sub) < 20 or sub["g"].nunique() < 2:
                    rows_t3.append({"pipe": pipe_name, "gene": gene, "model": mk,
                                    "rho": np.nan, "p": np.nan, "n": int(len(sub))})
                    continue
                y = ranknorm(sub["b"].values)
                x = ranknorm(sub["g"].values)
                C = (np.column_stack([np.ones(len(sub))] +
                                     [ranknorm(sub[c].values) for c in COVMAP[mk]])
                     if COVMAP[mk] else np.ones((len(sub), 1)))
                rr, pp = partial_pearson(y, x, C)
                rows_t3.append({"pipe": pipe_name, "gene": gene, "model": mk,
                                "rho": float(rr), "p": float(pp), "n": int(len(sub))})
    log(f"  T3 行数 {len(rows_t3)}（Sm_core：GBM 期望负相关，LGG 期望不显著）")


    # ---- 事件数下采样对照：把 SplAdder 高变异事件抽到 SpliceSeq 的同等规模 ----
    # 目的：silhouette / BIC / PAC 都对「样本数 : 特征数」比敏感，而两套管道的
    # 高变异事件数相差 3-7 倍。若抽稀后 SplAdder 也呈现 SpliceSeq 的指标，
    # 则差异应归因于维度而非球磨不生物学。
    n_target = min(Xs_.shape[1], Xa_.shape[1])
    rng = np.random.default_rng(SEED)
    ds_sil, ds_bicpref, ds_pac = [], [], []
    N_REP = 20
    log(f"  下采样对照：SplAdder {Xa_.shape[1]} 事件 -> {n_target}，重复 {N_REP} 次")
    for rep in range(N_REP):
        cols = rng.choice(Xa_.shape[1], size=n_target, replace=False)
        Xd = Xa_[:, cols]
        lab = KMeans(n_clusters=2, n_init=50, random_state=SEED).fit_predict(Xd)
        ds_sil.append(float(silhouette_score(Xd, lab)))
        b1 = GaussianMixture(1, covariance_type="full", random_state=SEED).fit(Xd).bic(Xd)
        b2 = GaussianMixture(2, covariance_type="full", random_state=SEED).fit(Xd).bic(Xd)
        ds_bicpref.append("k1" if b1 < b2 else "k2")
        if rep < 3:
            ds_pac.append(pac_monti(Xd, k=2, n_iter=60))
    ds = {"n_target": int(n_target),
          "sil_mean": float(np.mean(ds_sil)), "sil_sd": float(np.std(ds_sil)),
          "bic_k2_frac": float(np.mean([b == "k2" for b in ds_bicpref])),
          "bic_k1_frac": float(np.mean([b == "k1" for b in ds_bicpref])),
          "pac_mean": float(np.mean(ds_pac)) if ds_pac else np.nan}
    log(f"    下采样结果：sil={ds['sil_mean']:.3f}+/-{ds['sil_sd']:.3f}, "
        f"BIC偏好k2比例={ds['bic_k2_frac']:.2f}, PAC={ds['pac_mean']:.3f}")

    return {"cancer": cancer, "n_common": len(common),
            "n_event_SplAdder": int(nA), "n_event_SpliceSeq": int(nS),
            "T1_rho_spearman": float(rho_t1_s), "T1_p_spearman": float(p_t1_s),
            "T1_rho_pearson": float(rho_t1_p), "T1_p_pearson": float(p_t1_p),
            "SA_sil": ma["sil"], "SA_pac": ma["pac"], "SA_bic": ma["bic_prefer"],
            "SA_pc1var": pa["pc1var"], "SA_pc1b": pa["pc1_burden"],
            "SS_sil": ms["sil"], "SS_pac": ms["pac"], "SS_bic": ms["bic_prefer"],
            "SS_pc1var": ps["pc1var"], "SS_pc1b": ps["pc1_burden"],
            "rho_PC1_cross": float(rho_pc1_pipe),
            "ds_target": int(n_target), "ds_sil_mean": ds["sil_mean"], "ds_sil_sd": ds["sil_sd"],
            "ds_bic_k1_frac": ds["bic_k1_frac"], "ds_pac_mean": ds["pac_mean"],
            "ds": ds,
            "t3": rows_t3}


def main():
    results = []
    for cancer in ("GBM", "LGG"):
        r = process(cancer)
        if r:
            results.append(r)
    if not results:
        log("无结果写出。")
        return

    # ---- 写 T1/T2 ----
    rowsT1 = []
    rowsT2 = []
    rowsT3 = []
    for r in results:
        c = r["cancer"]
        rowsT1.append({"分析项": "T1-1 burden 直接相关", "队列": c,
                    "SplAdder (主管道)": "—",
                    "SpliceSeq (第二管道)": f"Spearman ρ={r['T1_rho_spearman']:.3f} (p={r['T1_p_spearman']:.2e})",
                    "备注": f"n={r['n_common']}"})
        rowsT1.append({"分析项": "T1-2 PC1 跨管道相关", "队列": c,
                    "SplAdder (主管道)": "—",
                    "SpliceSeq (第二管道)": f"Spearman ρ={r['rho_PC1_cross']:.3f}",
                    "备注": f"n={r['n_common']}"})
        rowsT2.append({"分析项": "T2-1 k=2 轮廓系数(silhouette)", "队列": c,
                    "SplAdder (主管道)": round(r["SA_sil"], 3),
                    "SpliceSeq (第二管道)": round(r["SS_sil"], 3),
                    "备注": f"事件 {r['n_event_SplAdder']} vs {r['n_event_SpliceSeq']}"})
        rowsT2.append({"分析项": "T2-2 修正版 PAC", "队列": c,
                    "SplAdder (主管道)": round(r["SA_pac"], 3),
                    "SpliceSeq (第二管道)": round(r["SS_pac"], 3), "备注": ""})
        rowsT2.append({"分析项": "T2-3 GMM-BIC 偏好", "队列": c,
                    "SplAdder (主管道)": r["SA_bic"], "SpliceSeq (第二管道)": r["SS_bic"], "备注": ""})
        rowsT2.append({"分析项": "T2-4 PC1 解释率", "队列": c,
                    "SplAdder (主管道)": round(r["SA_pc1var"], 3),
                    "SpliceSeq (第二管道)": round(r["SS_pc1var"], 3), "备注": ""})
        rowsT2.append({"分析项": "T2-5 PC1~burden ρ", "队列": c,
                    "SplAdder (主管道)": round(r["SA_pc1b"], 3),
                    "SpliceSeq (第二管道)": round(r["SS_pc1b"], 3), "备注": ""})
        rowsT2.append({"分析项": "T2-6 下采样对照(sil)", "队列": c,
                    "SplAdder (主管道)": f"{r['ds_sil_mean']:.3f}±{r['ds_sil_sd']:.3f}",
                    "SpliceSeq (第二管道)": round(r["SS_sil"], 3),
                    "备注": f"抽到 {r['ds_target']} 事件"})
        rowsT2.append({"分析项": "T2-7 下采样对照(BIC 偏好 k1 比例)", "队列": c,
                    "SplAdder (主管道)": round(r["ds_bic_k1_frac"], 2),
                    "SpliceSeq (第二管道)": r["SS_bic"], "备注": ""})
        rowsT2.append({"分析项": "T2-8 下采样对照(PAC)", "队列": c,
                    "SplAdder (主管道)": round(r["ds_pac_mean"], 3),
                    "SpliceSeq (第二管道)": round(r["SS_pac"], 3), "备注": ""})
        rowsT2.append({"分析项": "T2-9 下采样对照(SplAdder 事件抽稀)", "队列": c,
                       "SplAdder (主管道)": f"sil={r['ds']['sil_mean']:.3f}+/-{r['ds']['sil_sd']:.3f}, BIC偏好k2比例={r['ds']['bic_k2_frac']:.2f}",
                       "SpliceSeq (第二管道)": f"sil={r['SS_sil']:.3f}, BIC={r['SS_bic']}",
                       "备注": f"SplAdder {r['n_event_SplAdder']}->{r['ds']['n_target']} 事件"})
        for t in r["t3"]:
            rowsT3.append({"队列": c, "panel": t["gene"], "model": t["model"],
                           "SplAdder_rho": t["rho"], "SplAdder_p": t["p"],
                           "SpliceSeq_rho": None, "SpliceSeq_p": None, "n": t["n"]})
    # T3 两条管道配对
    rowsT3_pair = []
    for c in sorted({r["cancer"] for r in results}):
        rr = [r for r in results if r["cancer"] == c][0]
        a = {(t["gene"], t["model"]): t for t in rr["t3"] if t["pipe"] == "SplAdder"}
        b = {(t["gene"], t["model"]): t for t in rr["t3"] if t["pipe"] == "SpliceSeq"}
        for k in sorted(a):
            if k in b:
                rowsT3_pair.append({"分析项": f"T3 {k[0]} ({k[1]})", "队列": c,
                                    "SplAdder (主管道)": f"ρ={a[k]['rho']:.3f} (p={a[k]['p']:.2e})",
                                    "SpliceSeq (第二管道)": f"ρ={b[k]['rho']:.3f} (p={b[k]['p']:.2e})",
                                    "备注": f"n={a[k]['n']}"})

    T1 = pd.DataFrame(rowsT1)
    T2 = pd.DataFrame(rowsT2)
    T3p = pd.DataFrame(rowsT3_pair)
    T1.to_csv(os.path.join(OUT, "T1_results.csv"), index=False)
    T2.to_csv(os.path.join(OUT, "T2_cross_pipeline_results.csv"), index=False)
    T3p.to_csv(os.path.join(OUT, "T3_sm_results.csv"), index=False)

    log(f"\n{'='*72}\n【T1】\n{T1.to_string(index=False)}")
    log(f"\n{'='*72}\n【T2】\n{T2.to_string(index=False)}")
    log(f"\n{'='*72}\n【T3】\n{T3p.to_string(index=False)}")
    log(f"\n{'='*72}\n输出目录：{OUT}")


if __name__ == "__main__":
    main()
