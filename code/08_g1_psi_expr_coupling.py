#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 — 在 TCGA-GBM 内部复现「内含子滞留 → mRNA 下调」（NMD 链条本地化）

目的：摆脱对 GSE227187（乳腺癌）的依赖，在胶质瘤自身数据里建立 PSI–表达耦合。

分析 A（连续口径）：每个 IR 事件，跨 151 样本算 ρ(PSI, 宿主基因表达)
    · 主指标：负相关事件比例（NMD 预期 > 50%）
    · 零分布：置换样本标签（打乱 PSI 侧），保留表达矩阵的协方差结构
    · 混杂控制：偏相关，剔除全基因平均表达（MeanExprAll）的影响

分析 B（组间口径，与 GSE227187 逻辑一致）：
    RI-high vs RI-low 分组 → 每事件 ΔPSI、每基因 Δ表达
    · 检验：滞留升高的事件，其宿主基因表达是否更易下调
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
BASE = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = f"{BASE}/output"
FIG = f"{BASE}/figures"
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

MIN_DETECT = 0.90
N_PERM = 2000
SEED = 20260917


def barcode(c):
    return str(c).strip().strip('"')[:15]


def rank_rows(M):
    """按行做秩转换（NaN 保持 NaN）"""
    return M.apply(lambda r: r.rank(), axis=1)


def row_corr(A, B):
    """逐行 Pearson（A、B 形状相同，已中心化前会在此处理）"""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    a = A - np.nanmean(A, axis=1, keepdims=True)
    b = B - np.nanmean(B, axis=1, keepdims=True)
    a = np.nan_to_num(a)
    b = np.nan_to_num(b)
    num = (a * b).sum(axis=1)
    den = np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))
    return num / np.where(den == 0, np.nan, den)


def partial_out(R, cov):
    """把每行对协变量 cov（一维）做线性回归取残差"""
    c = cov - cov.mean()
    denom = (c * c).sum()
    R = np.asarray(R, dtype=float)
    Rc = R - np.nanmean(R, axis=1, keepdims=True)
    Rc = np.nan_to_num(Rc)
    beta = (Rc @ c) / denom
    return Rc - np.outer(beta, c)


def main():
    print("=" * 74)
    print("G1  TCGA-GBM 内部 PSI–表达耦合（NMD 链条本地化）")
    print("=" * 74)

    # ---------- 1. PSI ----------
    psi = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t",
                      index_col=0, compression="gzip", low_memory=False)
    meta = psi.iloc[:, :14].copy()
    vals = psi.iloc[:, 14:].copy()
    vals.columns = [barcode(c) for c in vals.columns]
    tumor = [c for c in vals.columns if c[13:15] == "01"]
    vals = vals[tumor].loc[:, ~vals[tumor].columns.duplicated()]
    P = vals.apply(pd.to_numeric, errors="coerce")
    print(f"[1] PSI: {P.shape[0]:,} 事件 × {P.shape[1]} 样本列")

    # ---------- 2. 表达 ----------
    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0,
                       compression="gzip")
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[~expr.index.duplicated()]
    common = [c for c in P.columns if c in set(expr.columns)]
    P = P[common]
    E = expr[common]
    print(f"[2] 与表达谱交集样本: {len(common)} 例；表达基因 {E.shape[0]:,}")
    print(f"    表达值范围: {np.nanmin(E.values):.2f} – {np.nanmax(E.values):.2f}")

    # ---------- 3. 事件过滤 ----------
    is_ir = (meta["SpliceType"] == "IR").values
    detect = P.notna().sum(axis=1).values / P.shape[1]
    Psd = P.std(axis=1, skipna=True).values
    gene = meta["GeneSymbol"].astype(str).values
    in_expr = np.array([g in set(E.index) for g in gene])
    keep = is_ir & (detect >= MIN_DETECT) & in_expr & (Psd > 0)
    print(f"[3] IR 事件 {int(is_ir.sum()):,} → 可分析（检出≥{MIN_DETECT}、基因在表达谱中、有变异）: {int(keep.sum()):,}")

    Pk = P.loc[keep, :]
    gk = pd.Index(gene[keep])
    Ek = E.loc[gk, :]                     # 与事件逐行对齐的宿主基因表达
    print(f"    涉及宿主基因: {gk.nunique():,}")

    # 缺失值：用事件中位数填补（秩转换后影响极小），并单独做敏感性分析
    P_fill = Pk.T.fillna(Pk.median(axis=1)).T
    complete = (Pk.notna().sum(axis=1) == Pk.shape[1]).values
    print(f"    完全无缺失的事件: {int(complete.sum()):,}（用于敏感性分析）")

    # ---------- 4. 分析 A：跨样本 Spearman ----------
    print("\n[4] 分析 A：跨样本 ρ(PSI, 宿主基因表达)")
    Rp = rank_rows(P_fill).values
    Re = rank_rows(Ek).values
    rho = row_corr(Rp, Re)

    n_neg = int(np.sum(rho < 0))
    n_ok = int(np.sum(~np.isnan(rho)))
    frac_neg = n_neg / n_ok
    bt = stats.binomtest(n_neg, n_ok, 0.5, alternative="two-sided")
    w = stats.wilcoxon(rho[~np.isnan(rho)], alternative="two-sided")
    print(f"    可算事件 {n_ok:,}｜负相关 {n_neg:,}（{frac_neg:.1%}）")
    print(f"    二项检验 vs 50%: p = {bt.pvalue:.3g}")
    print(f"    中位 ρ = {np.nanmedian(rho):+.4f}｜Wilcoxon（分布是否偏离 0）p = {w.pvalue:.3g}")

    # 置换零分布：打乱 PSI 侧的样本标签
    rng = np.random.default_rng(SEED)
    n_s = Rp.shape[1]
    null_frac, null_med = [], []
    Rpc = Rp - Rp.mean(axis=1, keepdims=True)
    Rec = Re - Re.mean(axis=1, keepdims=True)
    for _ in range(N_PERM):
        perm = rng.permutation(n_s)
        Rp2 = Rpc[:, perm]
        r2 = row_corr(Rp2, Rec)
        null_frac.append(np.nanmean(r2 < 0))
        null_med.append(np.nanmedian(r2))
    null_frac = np.array(null_frac)
    null_med = np.array(null_med)
    p_perm = (np.sum(np.abs(null_frac - 0.5) >= abs(frac_neg - 0.5)) + 1) / (N_PERM + 1)
    p_perm_med = (np.sum(np.abs(null_med) >= abs(np.nanmedian(rho))) + 1) / (N_PERM + 1)
    print(f"    置换零分布（{N_PERM} 次）：负相关比例 {null_frac.mean():.4f} ± {null_frac.std():.4f}"
          f"｜中位 ρ {null_med.mean():+.4f} ± {null_med.std():.4f}")
    print(f"    置换 p（负相关比例）= {p_perm:.4g}｜置换 p（中位 ρ）= {p_perm_med:.4g}")

    # ---------- 5. 混杂控制：剔除全基因平均表达 ----------
    print("\n[5] 偏相关：控制全基因平均表达（MeanExprAll）")
    mean_expr = E.mean(axis=0).values           # 每样本的平均表达
    mean_rank = stats.rankdata(mean_expr)
    mean_rank = (mean_rank - mean_rank.mean())
    Rp_res = partial_out(Rp, mean_rank)
    Re_res = partial_out(Re, mean_rank)
    rho_par = row_corr(Rp_res, Re_res)
    n_neg_p = int(np.sum(rho_par < 0))
    frac_neg_p = n_neg_p / np.sum(~np.isnan(rho_par))
    bt_p = stats.binomtest(n_neg_p, int(np.sum(~np.isnan(rho_par))), 0.5, alternative="two-sided")
    print(f"    负相关比例 {frac_neg_p:.1%}（{n_neg_p:,}）｜二项 p = {bt_p.pvalue:.3g}"
          f"｜中位 ρ = {np.nanmedian(rho_par):+.4f}")

    # ---------- 6. 敏感性：完全无缺失事件 ----------
    rho_c = rho[complete]
    fc = np.nanmean(rho_c < 0)
    print(f"\n[6] 敏感性（仅完全无缺失事件 n={int(complete.sum()):,}）："
          f"负相关 {fc:.1%}，中位 ρ = {np.nanmedian(rho_c):+.4f}")

    # ---------- 7. 分析 B：RI-high vs RI-low 组间 ----------
    print("\n[7] 分析 B：RI-high vs RI-low 组间 ΔPSI 与 Δ表达")
    ann = pd.read_csv(f"{OUT}/sample_RI_group_annotated.tsv", sep="\t", index_col=0)
    ann.index = [str(i)[:15] for i in ann.index]
    hi = [c for c in Pk.columns if ann.loc[c, "RI_group"] == "RI-high"]
    lo = [c for c in Pk.columns if ann.loc[c, "RI_group"] == "RI-low"]
    print(f"    RI-high {len(hi)} 例 / RI-low {len(lo)} 例")

    dPSI = Pk[hi].mean(axis=1) - Pk[lo].mean(axis=1)
    dExpr = Ek[hi].mean(axis=1) - Ek[lo].mean(axis=1)

    up = (dPSI > 0).values
    dn = (dPSI < 0).values
    mw = stats.mannwhitneyu(dExpr[up], dExpr[dn], alternative="two-sided")
    print(f"    滞留升高事件 {int(up.sum()):,}｜降低 {int(dn.sum()):,}")
    print(f"    Δ表达中位数：滞留升高组 {np.median(dExpr[up]):+.4f} vs 降低组 {np.median(dExpr[dn]):+.4f}")
    print(f"    Mann-Whitney p = {mw.pvalue:.3g}")
    rho_d = stats.spearmanr(dPSI.values, dExpr.values)
    print(f"    ρ(ΔPSI, Δ表达) = {rho_d.statistic:+.4f}, p = {rho_d.pvalue:.3g}")

    # 组间置换零分布
    rng2 = np.random.default_rng(SEED + 1)
    labs = np.array([1] * len(hi) + [0] * len(lo))
    allc = list(hi) + list(lo)
    Psub = Pk[allc].values
    Esub = Ek[allc].values
    null_mw, null_rd = [], []
    for _ in range(N_PERM):
        pm = rng2.permutation(len(labs))
        d1 = Psub[:, pm][:, labs == 1].mean(axis=1) - Psub[:, pm][:, labs == 0].mean(axis=1)
        d2 = Esub[:, pm][:, labs == 1].mean(axis=1) - Esub[:, pm][:, labs == 0].mean(axis=1)
        u = d1 > 0
        if u.sum() > 10 and (~u).sum() > 10:
            null_mw.append(np.median(d2[u]) - np.median(d2[~u]))
            null_rd.append(stats.spearmanr(d1, d2).statistic)
    null_mw = np.array(null_mw)
    null_rd = np.array(null_rd)
    obs_mw = np.median(dExpr[up]) - np.median(dExpr[dn])
    p_mw = (np.sum(np.abs(null_mw) >= abs(obs_mw)) + 1) / (len(null_mw) + 1)
    p_rd = (np.sum(np.abs(null_rd) >= abs(rho_d.statistic)) + 1) / (len(null_rd) + 1)
    print(f"    置换零分布：Δ表达差 {np.nanmean(null_mw):+.4f} ± {np.nanstd(null_mw):.4f}"
          f" → p = {p_mw:.4g}")
    print(f"                 ρ(ΔPSI,Δ表达) {np.nanmean(null_rd):+.4f} ± {np.nanstd(null_rd):.4f}"
          f" → p = {p_rd:.4g}")

    # ---------- 8. 按变异度分层 ----------
    print("\n[8] 按 PSI 变异度分层（高变异事件是否信号更强）")
    sd = Pk.std(axis=1).values
    qs = np.quantile(sd, [0.25, 0.5, 0.75])
    for lab, m in [("Q1 低变异", sd <= qs[0]), ("Q2–Q3", (sd > qs[0]) & (sd <= qs[2])),
                   ("Q4 高变异", sd > qs[2])]:
        sub = rho[m]
        sub = sub[~np.isnan(sub)]
        print(f"    {lab:12s} n={len(sub):6,}  负相关 {np.mean(sub < 0):.1%}  中位 ρ = {np.median(sub):+.4f}")

    # ---------- 9. 输出 ----------
    res = pd.DataFrame({
        "event": Pk.index, "gene": gk.values,
        "stdPSI": sd, "meanPSI": Pk.mean(axis=1).values,
        "rho_psi_expr": rho, "rho_partial": rho_par,
        "dPSI": dPSI.values, "dExpr": dExpr.values,
    })
    res.to_csv(f"{OUT}/G1_psi_expr_coupling.tsv", sep="\t", index=False)
    print(f"\n[9] 已输出 {OUT}/G1_psi_expr_coupling.tsv")

    summ = pd.DataFrame([
        ["A. 可分析 IR 事件", n_ok],
        ["A. 负相关事件比例", f"{frac_neg:.1%}"],
        ["A. 中位 rho", f"{np.nanmedian(rho):+.4f}"],
        ["A. 二项检验 p", f"{bt.pvalue:.3g}"],
        ["A. 置换 p（比例）", f"{p_perm:.4g}"],
        ["A. 置换 p（中位 rho）", f"{p_perm_med:.4g}"],
        ["A. 偏相关负相关比例", f"{frac_neg_p:.1%}"],
        ["A. 偏相关二项 p", f"{bt_p.pvalue:.3g}"],
        ["A. 敏感性（无缺失）负相关比例", f"{fc:.1%}"],
        ["B. 滞留升高事件数", int(up.sum())],
        ["B. 滞留降低事件数", int(dn.sum())],
        ["B. Δ表达中位数（升高组）", f"{np.median(dExpr[up]):+.4f}"],
        ["B. Δ表达中位数（降低组）", f"{np.median(dExpr[dn]):+.4f}"],
        ["B. Mann-Whitney p", f"{mw.pvalue:.3g}"],
        ["B. 置换 p", f"{p_mw:.4g}"],
        ["B. rho(ΔPSI, Δ表达)", f"{rho_d.statistic:+.4f}"],
        ["B. rho 置换 p", f"{p_rd:.4g}"],
    ], columns=["指标", "值"])
    summ.to_csv(f"{OUT}/G1_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
