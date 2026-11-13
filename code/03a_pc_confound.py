#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3a — RI 轴的混淆因素核查（最关键的一步）
目的：回答「PC1 / RI 负荷到底代表什么？」
若 RI 轴只是增殖、测序深度或肿瘤纯度的替身，整个设计作废。必须逐项排查。
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats as st

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
SEED = 20260917


def barcode(c):
    return str(c).strip().strip('"')[:15]


PANELS = {
    # 增殖（用于排除「RI 轴 = 增殖」）
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "BUB1",
                      "BUB1B", "AURKA", "AURKB", "PLK1", "CENPA", "CENPF", "RRM2",
                      "TYMS", "MCM2", "MCM6", "E2F1", "FOXM1", "TTK"],
    # 肿瘤纯度代理（胶质瘤高纯度 → 肿瘤细胞基因表达高）
    "TumorPurity_GBM_marker": ["GFAP", "SOX2", "OLIG2", "NES", "EGFR", "CD44",
                               "VIM", "PTPRZ1", "BCAN", "SPARC"],
    # 正常脑/神经元成分（纯度反向指标）
    "Neuronal_normal_brain": ["RBFOX3", "SYT1", "SNAP25", "GRIN1", "GABRA1",
                              "MAP2", "NEFL", "NEFM", "STMN2", "SYP"],
    # 剪接体/U1 snRNP 组分（验证 RI 轴的机制来源）
    "U1_snRNP": ["SNRPC", "SNRPA", "SNRPB", "SNRPD1", "SNRPD2", "SNRPD3",
                 "SNRPE", "SNRPF", "SNRPG", "SNRNP70", "SNRNP200", "PRPF8"],
    "Spliceosome_core": ["SF3B1", "SF3B2", "SF3A1", "U2AF1", "U2AF2", "SRSF1",
                         "SRSF2", "SRSF3", "HNRNPA1", "HNRNPC", "PTBP1", "RBM39"],
    # 干扰素/免疫（用于后续比较）
    "IFN_gamma_Ayers": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG",
                        "ISG20", "GBP1", "CD274", "CXCR6"],
    "Cytolytic": ["GZMA", "PRF1"],
    "Stromal": ["FAP", "POSTN", "COL1A1", "COL1A2", "COL3A1", "FN1", "DCN",
                "LUM", "THBS2", "ACTA2", "PDGFRB", "TIMP1"],
    "GliomaStemCell": ["SOX2", "NANOG", "POU5F1", "PROM1", "NES", "ALDH1A3",
                       "OLIG2", "SALL2", "BMI1", "KLF4", "MYC", "FUT4",
                       "ITGA6", "THY1", "PTPRZ1"],
}


def zscore_rows(df):
    """按基因（行）做 z 分数，返回基因 × 样本"""
    mu = df.mean(axis=1)
    sd = df.std(axis=1)
    return df.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)


def main():
    print("=" * 72)
    print("Step 3a  RI 轴混淆因素核查")
    print("=" * 72)
    d = np.load(f"{OUT}/output/cluster_input.npz", allow_pickle=True)
    M = d["M"]
    samples = [str(s) for s in d["samples"]]
    samp = pd.read_csv(f"{OUT}/output/sample_RI_group.tsv", sep="\t", index_col=0)

    # ---------- 1. PC1–PC5 ----------
    Z = (M - M.mean(0)) / M.std(0)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=5, random_state=SEED)
    P = pca.fit_transform(Z)
    for i in range(5):
        samp[f"RI_PC{i+1}"] = P[:, i]
    # PC1 方向对齐：使其与 RI 负荷正相关（便于解释）
    r0 = st.spearmanr(samp["RI_PC1"], samp["RI_burden"]).statistic
    if r0 < 0:
        for i in range(5):
            samp[f"RI_PC{i+1}"] *= -1
    print(f"[1] PC1 解释方差 {pca.explained_variance_ratio_[0]*100:.1f}%; "
          f"与 RI 负荷 Spearman rho = "
          f"{st.spearmanr(samp['RI_PC1'], samp['RI_burden']).statistic:.3f}")

    # ---------- 2. 表达谱模块分数 ----------
    expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0,
                       compression="gzip")
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[:, ~expr.columns.duplicated()]
    expr = expr[[c for c in samples if c in expr.columns]]
    print(f"[2] 表达谱: {expr.shape[0]:,} 基因 × {expr.shape[1]} 样本")

    expr_g = expr.groupby(level=0).mean() if expr.index.duplicated().any() else expr
    Zg = zscore_rows(expr_g)

    scores = pd.DataFrame(index=expr_g.columns)
    used = {}
    for name, genes in PANELS.items():
        have = [g for g in genes if g in Zg.index]
        used[name] = have
        if len(have) >= 3:
            scores[name] = Zg.loc[have].mean(axis=0)
        print(f"    {name:<26} {len(have)}/{len(genes)} 基因")
    # 文库复杂度代理：所有基因平均表达（≈ log 总表达量）
    scores["MeanExprAll"] = expr_g.mean(axis=0)
    # 线粒体基因比例（组织/RNA 质量代理）
    mt = [g for g in expr_g.index if str(g).startswith("MT-")]
    scores["MitoFrac"] = expr_g.loc[mt].mean(axis=0) if mt else np.nan

    # ---------- 3. 与 RI 负荷 / PC1 的相关 ----------
    tab = samp.join(scores, how="left")
    targets = ["RI_burden", "RI_PC1", "RI_PC2"]
    print("\n[3] 混淆核查（Spearman，n={}）:".format(len(tab)))
    print(f"    {'变量':<26} {'vs RI_burden':>14} {'p':>10} {'vs PC1':>10} {'p':>10}")
    rows = []
    for col in list(PANELS.keys()) + ["MeanExprAll", "MitoFrac"]:
        if col not in tab.columns:
            continue
        v = pd.to_numeric(tab[col], errors="coerce")
        ok = v.notna() & tab["RI_burden"].notna()
        if ok.sum() < 20:
            continue
        r1, p1 = st.spearmanr(tab.loc[ok, "RI_burden"], v[ok])
        r2, p2 = st.spearmanr(tab.loc[ok, "RI_PC1"], v[ok])
        print(f"    {col:<26} {r1:>14.3f} {p1:>10.2e} {r2:>10.3f} {p2:>10.2e}")
        rows.append({"variable": col, "rho_burden": r1, "p_burden": p1,
                     "rho_PC1": r2, "p_PC1": p2, "n": int(ok.sum()),
                     "n_genes_used": len(used.get(col, []))})
        # 若与纯度/增殖高度相关，标记
    corr = pd.DataFrame(rows)
    corr.to_csv(f"{OUT}/output/step3a_confound_corr.tsv", sep="\t", index=False)

    flag = corr[(corr.variable.isin(["Proliferation", "Neuronal_normal_brain",
                                     "TumorPurity_GBM_marker", "MeanExprAll",
                                     "MitoFrac"])) &
                (corr.rho_burden.abs() > 0.5)]
    if len(flag):
        print("\n    [!] 以下混淆变量与 RI 负荷 |rho|>0.5，必须在论文中校正/说明：")
        print(flag[["variable", "rho_burden", "p_burden"]].to_string(index=False))
    else:
        print("\n    [OK] 增殖/纯度/表达量代理与 RI 负荷均 |rho|<0.5")

    # ---------- 4. 偏相关：剔除增殖后 RI 是否仍与免疫相关 ----------
    print("\n[4] 偏相关（控制增殖与纯度）:")
    def partial_spearman(x, y, ctrl_names):
        sub = tab[[x, y] + ctrl_names].apply(pd.to_numeric, errors="coerce").dropna()
        if len(sub) < 30:
            return np.nan, np.nan, 0
        rx = st.rankdata(sub[x]); ry = st.rankdata(sub[y])
        C = np.column_stack([np.ones(len(sub))] +
                            [st.rankdata(sub[c]) for c in ctrl_names])
        bx = np.linalg.lstsq(C, rx, rcond=None)[0]
        by = np.linalg.lstsq(C, ry, rcond=None)[0]
        return st.pearsonr(rx - C @ bx, ry - C @ by) + (len(sub),)
    for y in ["IFN_gamma_Ayers", "Cytolytic", "GliomaStemCell", "Stromal"]:
        if y not in tab.columns:
            continue
        r, p, n = partial_spearman("RI_burden", y, ["Proliferation",
                                                    "Neuronal_normal_brain"])
        print(f"    RI_burden ~ {y:<18} 偏 rho={r:+.3f}, p={p:.2e}, n={n}")

    # ---------- 5. 分组差异（除负荷外还有何不同）----------
    if "RI_PC1" in tab.columns:
        from scipy.stats import mannwhitneyu
        print("\n[5] RI-high vs RI-low 的关键分数差异:")
        for col in ["RI_PC1", "Proliferation", "GliomaStemCell", "IFN_gamma_Ayers",
                    "Cytolytic", "Stromal", "Neuronal_normal_brain", "MeanExprAll"]:
            if col not in tab.columns:
                continue
            a = pd.to_numeric(tab.loc[tab.RI_group == "RI-high", col], errors="coerce").dropna()
            b = pd.to_numeric(tab.loc[tab.RI_group == "RI-low", col], errors="coerce").dropna()
            if len(a) < 5 or len(b) < 5:
                continue
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            print(f"    {col:<26} high={a.mean():+7.3f} low={b.mean():+7.3f} p={p:.3e}")

    tab.to_csv(f"{OUT}/output/sample_RI_group_annotated.tsv", sep="\t")
    with open(f"{OUT}/output/step3a_panels_used.json", "w") as f:
        json.dump(used, f, ensure_ascii=False, indent=2)
    print("\n完成 -> sample_RI_group_annotated.tsv / step3a_confound_corr.tsv")


if __name__ == "__main__":
    main()
