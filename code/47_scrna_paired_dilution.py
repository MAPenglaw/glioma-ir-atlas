#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 47 — 单细胞样本内配对 + "浸润稀释"假说定量分解（方案 C 第一部分）

科学问题
--------
bulk 层面发现 IR burden 与免疫/基质浸润呈负相关（LGG 17/17，GBM 髓系为主）。
一个必须排除的备择解释是 **浸润稀释（infiltration dilution）**：

    若免疫细胞自身的 IR 轴分数低于肿瘤细胞，则免疫浸润比例 f 越高，
    bulk 观测值 = (1-f)·T + f·I 就越低 —— 负相关完全由**组成**产生，
    与肿瘤细胞内在生物学无关。

本脚本在 44 例单细胞样本内做配对，把这个问题拆成可检验的三段：

  T_s  = 样本 s 中恶性细胞（Tumor）的平均 IR 轴代理分
  I_s  = 样本 s 中免疫细胞（Myeloid/T_cell/Immune_other 按实际比例加权）的平均分
  f_s  = 样本 s 中免疫细胞占比

  (1) 样本内配对 T_s − I_s 是否 ≠ 0          → 稀释在数学上是否成立
  (2) corr(T_s, f_s) 是否 ≠ 0                → 肿瘤**内在**关联是否存在
  (3) 反事实分解：固定 T 只让 f 变 / 固定 I 只让 f 变
                                              → 稀释与内在各贡献多少

深度混杂
--------
细胞类型间的 total_counts 差异会让"平均表达"失真。Fig6 现行口径只看了
RIaxis_pos（正相关基因集），而 pos 与 neg 在各细胞类型间同向变化，
提示存在细胞类型层面的整体表达偏移。因此本脚本同时报告：
  - 未校正（score = pos − neg，已含内部对照）
  - 深度校正残差（score ~ log10(total_counts)，全细胞 pooled）
  - 深度分层匹配（按 total_counts 十分位层内比较后加权）

输出（output/）
--------------
  sc_paired_cells.tsv      逐细胞（含 total_counts / 深度校正残差）
  sc_paired_sample.tsv     每样本 × 每区室的汇总（配对分析输入）
  sc_paired_stats.tsv      配对检验 + 稀释分解的全部统计量
"""

import gzip
import os
import re

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy import stats as st

SC = "/Users/peng/课题/省卫健委/sc_data/GSE182109"
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output"

IMMUNE = ["Myeloid", "T_cell", "Immune_other"]
MIN_CELLS = 20          # 每样本每个区室参与配对的最少细胞数
N_BOOT = 2000


# ---------------------------------------------------------------- 数据加载
def load_sample(prefix, sample_id, want):
    """读取一个 10x 样本：目标基因子集 + 每细胞 total_counts。"""
    feat = pd.read_csv(prefix + "_features.tsv.gz", sep="\t", header=None)
    syms = feat[1].astype(str).values
    pos = {}
    for i, s in enumerate(syms):
        if s in want and s not in pos:
            pos[s] = i
    if not pos:
        return None
    with gzip.open(prefix + "_matrix.mtx.gz", "rb") as f:
        M = sio.mmread(f).tocsr()
    totals = np.asarray(M.sum(axis=0)).ravel()
    idx = np.array(list(pos.values()))
    sub = M[idx, :].toarray().astype(np.float32)
    df = pd.DataFrame(sub.T, columns=list(pos.keys()))
    bc = pd.read_csv(prefix + "_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).values
    if len(bc) != len(df):
        return None
    df["total_counts"] = totals
    df["sample"] = sample_id
    return df


def annotate(norm, POS, NEG):
    """与 step 06 完全一致的细胞类型注释 + 轴分数。"""
    ptprc = norm["PTPRC"].values > 0
    tmark = (norm["CD3E"].values > 0) | (norm["CD8A"].values > 0)
    mmark = ((norm["AIF1"].values > 0) | (norm["CD68"].values > 0)
             | (norm["C1QA"].values > 0))
    tum = (norm["SOX2"].values > 0) | (norm["EGFR"].values > 0)
    norm["cell_type"] = np.where(
        ptprc & tmark, "T_cell",
        np.where(ptprc & mmark, "Myeloid",
        np.where(ptprc, "Immune_other",
        np.where(tum, "Tumor",
        np.where(norm["CLDN5"].values > 0, "Endothelial",
        np.where(norm["COL1A1"].values > 0, "Fibroblast",
        np.where(norm["PLP1"].values > 0, "Oligo_like", "Other")))))))

    def mscore(cols):
        cs = [c for c in cols if c in norm.columns]
        return norm[cs].mean(axis=1).values if cs else np.full(len(norm), np.nan)

    norm["RIaxis_pos"] = mscore(POS)
    norm["RIaxis_neg"] = mscore(NEG)
    norm["RIaxis_score"] = norm["RIaxis_pos"] - norm["RIaxis_neg"]
    return norm


def lineage(sid):
    if sid.startswith("LGG"):
        return "LGG"
    if sid.startswith("ndGBM"):
        return "ndGBM"
    if sid.startswith("rGBM"):
        return "rGBM"
    return "other"


def main():
    sig = pd.read_csv(os.path.join(OUT, "RI_axis_signature_genes.tsv"), sep="\t")
    POS = sig.head(100)["gene"].tolist()
    NEG = sig.tail(100)["gene"].tolist()
    MARKERS = ["PTPRC", "CD3E", "CD8A", "AIF1", "CD68", "C1QA", "SOX2", "EGFR",
               "CLDN5", "COL1A1", "PLP1", "SNRPC", "UQCRC1"]
    WANT = set(POS + NEG + MARKERS)

    prefixes = sorted({f[:-len("_barcodes.tsv.gz")] for f in os.listdir(SC)
                       if f.endswith("_barcodes.tsv.gz")})
    print(f"样本数: {len(prefixes)}", flush=True)

    dfs = []
    for n, p in enumerate(prefixes):
        sid = os.path.basename(p).split("_", 1)[1]
        try:
            d = load_sample(os.path.join(SC, p), sid, WANT)
            if d is not None:
                dfs.append(d)
        except Exception as e:
            print(f"  [跳过] {sid}: {e}", flush=True)
        if (n + 1) % 10 == 0:
            print(f"  [{n+1}/{len(prefixes)}]", flush=True)

    ad = pd.concat(dfs, axis=0)
    ad = ad[ad.total_counts > 500].copy()
    print("QC 后细胞数:", len(ad), flush=True)

    genes = [c for c in ad.columns if c not in ("total_counts", "sample")]
    norm = np.log1p(ad[genes].div(ad["total_counts"], axis=0) * 1e4)
    norm["total_counts"] = ad["total_counts"].values
    norm["sample"] = ad["sample"].values
    norm = annotate(norm, POS, NEG)

    norm["lineage"] = norm["sample"].map(lineage)
    # 病人 ID：去掉末尾的区域编号（-1/-A/-B ...）
    norm["patient"] = norm["sample"].str.replace(r"-[A-Z0-9]+$", "", regex=True)
    norm["log_umi"] = np.log10(norm["total_counts"].values)

    # ---- 深度诊断：各细胞类型的 UMI 分布 ----
    print("\n=== 深度诊断：各细胞类型 total_counts 中位数 ===")
    dep = norm.groupby("cell_type")["total_counts"].agg(
        ["median", "mean", "count"]).sort_values("median", ascending=False)
    print(dep.round(0).to_string())
    print("\nTumor vs 免疫 UMI Mann-Whitney p = %.2e" % st.mannwhitneyu(
        norm.loc[norm.cell_type == "Tumor", "total_counts"],
        norm.loc[norm.cell_type.isin(IMMUNE), "total_counts"]).pvalue)

    # ---- 深度校正残差：全细胞 pooled 回归 score ~ log10(UMI) ----
    X = np.column_stack([np.ones(len(norm)), norm["log_umi"].values])
    b = np.linalg.lstsq(X, norm["RIaxis_score"].values, rcond=None)[0]
    norm["score_adj"] = norm["RIaxis_score"].values - X @ b
    print(f"深度校正: score = {b[0]:+.4f} {b[1]:+.4f} * log10(UMI)")

    # ---- 深度分层匹配：按总 UMI 十分位分层 ----
    norm["depth_bin"] = pd.qcut(norm["total_counts"], 10, labels=False,
                                duplicates="drop")

    norm[["cell_type", "sample", "patient", "lineage", "RIaxis_pos",
          "RIaxis_neg", "RIaxis_score", "score_adj", "total_counts",
          "log_umi", "depth_bin", "SNRPC", "UQCRC1"]].to_csv(
        os.path.join(OUT, "sc_paired_cells.tsv"), sep="\t", index=False)
    print("\n-> sc_paired_cells.tsv")

    # ---------------------------------------------------- 每样本 × 区室汇总
    rows = []
    for (sid), g in norm.groupby("sample"):
        tum = g[g.cell_type == "Tumor"]
        imm = g[g.cell_type.isin(IMMUNE)]
        if len(tum) < MIN_CELLS or len(imm) < MIN_CELLS:
            continue
        # 免疫区室内部按实际比例加权（合成 bulk 时需要）
        ic = imm["cell_type"].value_counts()
        w = (ic / ic.sum())
        i_pos = float(np.average([imm.loc[imm.cell_type == k, "RIaxis_pos"].mean()
                                  for k in w.index], weights=w.values))
        i_neg = float(np.average([imm.loc[imm.cell_type == k, "RIaxis_neg"].mean()
                                  for k in w.index], weights=w.values))
        i_score = float(np.average([imm.loc[imm.cell_type == k, "RIaxis_score"].mean()
                                    for k in w.index], weights=w.values))
        i_adj = float(np.average([imm.loc[imm.cell_type == k, "score_adj"].mean()
                                  for k in w.index], weights=w.values))
        rows.append(dict(
            sample=sid, patient=g["patient"].iloc[0], lineage=g["lineage"].iloc[0],
            n_cells=len(g),
            n_tumor=len(tum), n_immune=len(imm),
            f_immune=len(imm) / len(g),
            T_pos=tum["RIaxis_pos"].mean(), T_neg=tum["RIaxis_neg"].mean(),
            T_score=tum["RIaxis_score"].mean(), T_adj=tum["score_adj"].mean(),
            I_pos=i_pos, I_neg=i_neg, I_score=i_score, I_adj=i_adj,
            umi_tumor=tum["total_counts"].median(),
            umi_immune=imm["total_counts"].median(),
        ))
    S = pd.DataFrame(rows)
    S["d_score"] = S["T_score"] - S["I_score"]        # 未校正（score 已含内部对照）
    S["d_adj"] = S["T_adj"] - S["I_adj"]              # 深度校正后
    S["d_pos"] = S["T_pos"] - S["I_pos"]              # Fig6 现行口径（无内部对照）
    S.to_csv(os.path.join(OUT, "sc_paired_sample.tsv"), sep="\t", index=False)
    print(f"-> sc_paired_sample.tsv  (n={len(S)} 样本可用于配对)")

    # ---------------------------------------------------------- 配对检验
    out = []

    def paired(col, label):
        v = S[col].dropna().values
        n = len(v)
        if n < 5:
            return
        w = st.wilcoxon(v, alternative="two-sided")
        t = st.ttest_1samp(v, 0)
        # 统一口径：报告**均值差**及其 bootstrap 95% CI（不再混用中位数）
        boot = np.random.default_rng(0).choice(v, (N_BOOT, n), replace=True).mean(axis=1)
        out.append(dict(metric=label, n=n,
                        mean_diff=float(v.mean()),
                        sd=float(v.std(ddof=1)),
                        median_diff=float(np.median(v)),
                        ci_lo=float(np.percentile(boot, 2.5)),
                        ci_hi=float(np.percentile(boot, 97.5)),
                        wilcoxon_p=float(w.pvalue),
                        t_p=float(t.pvalue),
                        pos_frac=float((v > 0).mean())))

    print("\n=== 样本内配对：恶性 vs 免疫（44 例）===")
    paired("d_score", "T-I score (internal control)")
    paired("d_adj", "T-I score (depth-adjusted)")
    paired("d_pos", "T-I pos (Fig6 current)")
    for r in out:
        print(f"  {r['metric']:<26} n={r['n']}  中位差={r['median_diff']:+.4f} "
              f"[{r['ci_lo']:+.4f},{r['ci_hi']:+.4f}]  Wilcoxon p={r['wilcoxon_p']:.3g} "
              f" T>0 占 {r['pos_frac']*100:.0f}%")

    # ------------------------------------------- 稀释假说：定量分解
    print("\n=== 稀释假说分解 ===")
    f = S["f_immune"].values
    T = S["T_score"].values
    I = S["I_score"].values

    bulk_obs = (1 - f) * T + f * I                    # 实际合成 bulk
    cf_dilute = (1 - f) * T.mean() + f * I            # 反事实：肿瘤内在固定
    cf_intrin = (1 - f) * T + f * I.mean()            # 反事实：免疫组成固定

    def sc(a, b):
        r, p = st.spearmanr(a, b)
        return float(r), float(p)

    rows2 = []
    for nm, v in [("合成 bulk（观测）", bulk_obs),
                  ("反事实：仅稀释（T 固定）", cf_dilute),
                  ("反事实：仅肿瘤内在（I 固定）", cf_intrin)]:
        r, p = sc(v, f)
        rows2.append(dict(quantity=nm, rho_vs_f=r, p=p))
        print(f"  {nm:<28} rho vs f = {r:+.3f} (p={p:.3g})")

    r_tf, p_tf = sc(T, f)
    r_if, p_if = sc(I, f)
    rows2.append(dict(quantity="肿瘤内在 T vs f", rho_vs_f=r_tf, p=p_tf))
    rows2.append(dict(quantity="免疫区室 I vs f", rho_vs_f=r_if, p=p_if))
    print(f"  {'肿瘤内在 T vs f':<28} rho = {r_tf:+.3f} (p={p_tf:.3g})")
    print(f"  {'免疫区室 I vs f':<28} rho = {r_if:+.3f} (p={p_if:.3g})")

    # 稀释能解释的比例：纯稀释相关 / 观测相关
    r_obs = rows2[0]["rho_vs_f"]
    r_dil = rows2[1]["rho_vs_f"]
    frac = (r_dil / r_obs) if abs(r_obs) > 1e-6 else np.nan
    rows2.append(dict(quantity="稀释可解释比例 (rho_dilute/rho_obs)",
                      rho_vs_f=frac, p=np.nan))
    print(f"\n  稀释可解释比例 = {frac:.3f}" if np.isfinite(frac)
          else "\n  稀释可解释比例: n/a")

    # ------------------------------- 方差分解 + 置换检验（稀释机制的产能上限）
    # 反事实 SD：纯稀释 / 纯肿瘤内在 各自能制造多少 bulk 变异
    f_v, T_v, I_v = S["f_immune"].values, S["T_score"].values, S["I_score"].values
    bulk_v = (1 - f_v) * T_v + f_v * I_v
    cf_d_v = (1 - f_v) * T_v.mean() + f_v * I_v
    cf_i_v = (1 - f_v) * T_v + f_v * I_v.mean()
    vdec = pd.DataFrame([
        dict(quantity="observed synthetic bulk", sd=float(bulk_v.std()), pct=1.0),
        dict(quantity="dilution only (T fixed)", sd=float(cf_d_v.std()),
             pct=float(cf_d_v.std() / bulk_v.std())),
        dict(quantity="tumour-intrinsic only (I fixed)", sd=float(cf_i_v.std()),
             pct=float(cf_i_v.std() / bulk_v.std())),
        dict(quantity="between-sample SD of T", sd=float(T_v.std()), pct=np.nan),
        dict(quantity="mean T-I difference", sd=float(S["d_score"].mean()), pct=np.nan),
    ])
    print("\n=== 方差分解（SD 幅度）===")
    print(vdec.round(4).to_string(index=False))

    # 置换：打乱 f 与 (T, I) 的配对 → T、I 都与 f 独立 = 纯稀释机制
    rng = np.random.default_rng(0)
    perm = np.empty(5000)
    for i in range(5000):
        fp = rng.permutation(f_v)
        perm[i] = st.spearmanr((1 - fp) * T_v + fp * I_v, fp).statistic
    # TCGA 实测参考值（45 号 M5，供对比）
    REF = {"LGG IFN-gamma": -0.384, "LGG CD8_T": -0.321, "GBM myeloid": -0.364}
    pr = {k: float(np.mean(perm <= v)) for k, v in REF.items()}
    print("\n=== 置换检验：纯稀释机制的 rho 分布（5000 次）===")
    print(f"  中位 {np.median(perm):+.3f}  95% 区间 "
          f"[{np.percentile(perm,2.5):+.3f}, {np.percentile(perm,97.5):+.3f}]")
    for k, v in pr.items():
        print(f"  达到 TCGA 实测 {k} ({v:+.3f}) 的概率 = {pr[k]*100:.2f}%")
    pd.DataFrame({"perm_rho": perm}).to_csv(os.path.join(OUT, "sc_dilution_perm.tsv"),
                                            sep="\t", index=False)
    vdec.to_csv(os.path.join(OUT, "sc_dilution_vardec.tsv"), sep="\t", index=False)
    pd.DataFrame([dict(reference=k, tcga_rho=v, p_dilution_ge=pr[k])
                  for k, v in REF.items()]).to_csv(
        os.path.join(OUT, "sc_dilution_permref.tsv"), sep="\t", index=False)

    # ------------------------------------------- 深度分层匹配（敏感性）
    print("\n=== 深度分层匹配（UMI 十分位层内比较）===")
    layer = []
    for b_, g in norm.groupby("depth_bin"):
        t_ = g.loc[g.cell_type == "Tumor", "RIaxis_score"]
        i_ = g.loc[g.cell_type.isin(IMMUNE), "RIaxis_score"]
        if len(t_) < 50 or len(i_) < 50:
            continue
        layer.append((b_, len(t_), len(i_), t_.mean() - i_.mean()))
    L = pd.DataFrame(layer, columns=["bin", "n_t", "n_i", "diff"])
    if len(L):
        wm = np.average(L["diff"], weights=L["n_t"] + L["n_i"])
        print(f"  层数={len(L)}  加权平均差={wm:+.4f}  "
              f"各层全负: {bool((L['diff']<0).all())}  各层全正: {bool((L['diff']>0).all())}")
        print(L.round(4).to_string(index=False))
        out.append(dict(metric="深度分层匹配 加权差", n=len(L), median_diff=float(wm),
                        ci_lo=np.nan, ci_hi=np.nan, wilcoxon_p=np.nan,
                        t_p=np.nan, pos_frac=float((L["diff"] > 0).mean())))
    else:
        print("  层内细胞不足，跳过")

    pd.DataFrame(out).to_csv(os.path.join(OUT, "sc_paired_stats.tsv"),
                             sep="\t", index=False)
    pd.DataFrame(rows2).to_csv(os.path.join(OUT, "sc_dilution_decomp.tsv"),
                               sep="\t", index=False)
    print("\n-> sc_paired_stats.tsv / sc_dilution_decomp.tsv")


if __name__ == "__main__":
    main()
