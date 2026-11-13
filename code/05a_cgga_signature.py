#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 5a — RI 轴的表达层面特征签名 + CGGA 外部验证（1,018 例）

逻辑：
  既然 TCGA-GBM 的 PSI 层面 RI 负荷没有预后价值，就需要问一个更有意义的问题：
  「RI 轴所对应的转录组状态，在更大的外部队列里是否有临床意义？」
  做法：在 TCGA-GBM 中找出表达水平上与 RI 负荷相关的基因（RI 轴表达签名），
        再把这个签名搬到 CGGA（693+325，含 WHO II–IV、IDH/1p19q/MGMT）里检验。
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats as st

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
CG = f"{OUT}/data/cgga"


def bc(c):
    return str(c).strip().strip('"')[:15]


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


# ---------- 1. TCGA-GBM: 求 RI 轴表达签名 ----------
print("=" * 72)
print("Step 5a  RI 轴表达签名 + CGGA 外部验证")
print("=" * 72)
tab = pd.read_csv(f"{OUT}/output/sample_RI_group_clinical.tsv", sep="\t", index_col=0)
burden = pd.to_numeric(tab["RI_burden"], errors="coerce")

expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0, compression="gzip")
expr.columns = [bc(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()]
expr = expr.groupby(level=0).mean()
common = [s for s in burden.index if s in expr.columns]
X = expr[common]
y = burden[common].values
print(f"[1] TCGA-GBM: {X.shape[0]:,} 基因 × {len(common)} 样本")

# 过滤低表达基因（50% 样本 > log2(1)=0）
keep = (X > 1).sum(axis=1) >= 0.5 * X.shape[1]
Xk = X[keep]
print(f"    通过表达过滤: {Xk.shape[0]:,} 基因")

Xv = Xk.values.astype(float)
ry = st.rankdata(y)
rx = np.apply_along_axis(st.rankdata, 1, Xv)
rx_c = rx - rx.mean(axis=1, keepdims=True)
ry_c = ry - ry.mean()
num = (rx_c * ry_c).sum(axis=1)
den = np.sqrt((rx_c ** 2).sum(axis=1) * (ry_c ** 2).sum())
rho = num / np.where(den == 0, np.nan, den)
n = len(y)
tstat = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
pval = 2 * st.t.sf(np.abs(tstat), df=n - 2)
qval = bh(pval)
sig = pd.DataFrame({"gene": Xk.index, "rho": rho, "p": pval, "q": qval})
sig = sig[sig.q < 0.05].sort_values("rho", ascending=False)
print(f"    FDR<0.05 的 RI 相关基因: {len(sig)} "
      f"(正相关 {int((sig.rho>0).sum())}, 负相关 {int((sig.rho<0).sum())})")
sig.to_csv(f"{OUT}/output/RI_axis_signature_genes.tsv", sep="\t", index=False)

NPOS, NNEG = 100, 100
pos = sig.head(NPOS)["gene"].tolist()
neg = sig.tail(NNEG)["gene"].tolist()
print(f"    签名取 rho 最高 {len(pos)} 个 + 最低 {len(neg)} 个基因")

# ---------- 2. 载入 CGGA ----------
def load_cgga(expr_f, clin_f, tag):
    cl = pd.read_csv(clin_f, sep="\t", index_col=0, dtype=str)
    cl.index = cl.index.str.strip()
    ex = pd.read_csv(expr_f, sep="\t", index_col=0)
    ex.index = [str(i).strip() for i in ex.index]
    ex = ex.groupby(level=0).mean()
    ex.columns = [c.strip() for c in ex.columns]
    ex = np.log2(ex + 1)
    cl = cl.loc[[i for i in cl.index if i in ex.columns]]
    ex = ex[cl.index]
    cl = cl.rename(columns=lambda c: c.strip())
    # 关键字段
    def g(*names):
        for nm in names:
            for c in cl.columns:
                if c.lower().startswith(nm.lower()):
                    return c
        return None
    out = pd.DataFrame(index=cl.index)
    out["Grade"] = cl[g("Grade")].str.extract(r"WHO\s*(I{1,3}V?|IV)", expand=False) \
        if g("Grade") else np.nan
    out["Grade_raw"] = cl[g("Grade")] if g("Grade") else np.nan
    out["Gender"] = cl[g("Gender")] if g("Gender") else np.nan
    out["Age"] = pd.to_numeric(cl[g("Age")], errors="coerce") if g("Age") else np.nan
    out["OS"] = pd.to_numeric(cl[g("OS")], errors="coerce") if g("OS") else np.nan
    cc = g("Censor")
    out["Censor"] = pd.to_numeric(cl[cc], errors="coerce") if cc else np.nan
    out["IDH"] = cl[g("IDH")] if g("IDH") else np.nan
    out["codel_1p19q"] = cl[g("1p19q")] if g("1p19q") else np.nan
    out["MGMT"] = cl[g("MGMT")] if g("MGMT") else np.nan
    out["Histology"] = cl[g("Histology")] if g("Histology") else np.nan
    out["cohort"] = tag
    # 签名打分
    zp = ((ex.loc[[x for x in pos if x in ex.index]] -
           ex.loc[[x for x in pos if x in ex.index]].mean(axis=1).values[:, None]) /
          ex.loc[[x for x in pos if x in ex.index]].std(axis=1).replace(0, np.nan).values[:, None])
    zn = ((ex.loc[[x for x in neg if x in ex.index]] -
           ex.loc[[x for x in neg if x in ex.index]].mean(axis=1).values[:, None]) /
          ex.loc[[x for x in neg if x in ex.index]].std(axis=1).replace(0, np.nan).values[:, None])
    out["RIaxis_score"] = zp.mean(axis=0).values - zn.mean(axis=0).values
    out["n_pos_used"] = len([x for x in pos if x in ex.index])
    out["n_neg_used"] = len([x for x in neg if x in ex.index])
    print(f"[2] {tag}: {ex.shape[0]:,} 基因 × {ex.shape[1]} 样本; "
          f"签名命中 {out['n_pos_used'].iloc[0]}+{out['n_neg_used'].iloc[0]} 基因; "
          f"OS 可用 {out['OS'].notna().sum()}")
    return out


c693 = load_cgga(f"{CG}/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
                 f"{CG}/CGGA.mRNAseq_693_clinical.20200506.txt", "CGGA_693")
c325 = load_cgga(f"{CG}/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
                 f"{CG}/CGGA_mRNAseq_325_clinical.20200506.txt"
                 if os.path.exists(f"{CG}/CGGA_mRNAseq_325_clinical.20200506.txt")
                 else f"{CG}/CGGA.mRNAseq_325_clinical.20200506.txt", "CGGA_325")

cg = pd.concat([c693, c325])
cg = cg[cg.OS.notna() & cg.Censor.notna()]
cg["event"] = cg["Censor"].astype(int)
cg["OS_months"] = cg["OS"] / 30.44
print(f"[3] 合并 CGGA: {len(cg)} 例（去除重复 ID: "
      f"{len(cg) - len(cg.index.unique())}）")
cg = cg[~cg.index.duplicated(keep="first")]
print(f"    最终: {len(cg)} 例 | 事件 {int(cg.event.sum())} | "
      f"IDH 已知 {cg.IDH.notna().sum()} | 分级已知 {cg.Grade.notna().sum()}")
print("    IDH 分布:", cg.IDH.value_counts(dropna=False).to_dict())
print("    分级分布:", cg.Grade.value_counts(dropna=False).to_dict())

cg.to_csv(f"{OUT}/output/cgga_riaxis_scores.tsv", sep="\t")

# ---------- 3. TCGA 内的一致性检查（签名 vs 负荷 / U1 / 干性）----------
sc = X.loc[[g for g in pos + neg if g in X.index]].copy()
zp = ((sc.loc[[g for g in pos if g in sc.index]] -
       sc.loc[[g for g in pos if g in sc.index]].mean(axis=1).values[:, None]) /
      sc.loc[[g for g in pos if g in sc.index]].std(axis=1).replace(0, np.nan).values[:, None])
zn = ((sc.loc[[g for g in neg if g in sc.index]] -
       sc.loc[[g for g in neg if g in sc.index]].mean(axis=1).values[:, None]) /
      sc.loc[[g for g in neg if g in sc.index]].std(axis=1).replace(0, np.nan).values[:, None])
tga_score = (zp.mean(axis=0) - zn.mean(axis=0))
r, p = st.spearmanr(tga_score.values, burden[tga_score.index].values)
print(f"\n[4] TCGA 内: 签名分数 vs RI 负荷 Spearman rho = {r:.3f} (p={p:.2e})")
gsva = pd.read_csv(f"{OUT}/output/gsva_scores.tsv", sep="\t", index_col=0)
for c in ["U1_snRNP", "Spliceosome_core", "Glioma_StemCell", "Proliferation",
          "Treg", "CD8_T_cell"]:
    if c in gsva.columns:
        rr, pp = st.spearmanr(tga_score.values, gsva.loc[tga_score.index, c].values)
        print(f"    签名 vs GSVA {c:<20} rho={rr:+.3f} (p={pp:.2e})")

with open(f"{OUT}/output/step5a_summary.json", "w") as f:
    json.dump({"n_sig_genes": int(len(sig)), "n_pos": len(pos), "n_neg": len(neg),
               "tcga_score_vs_burden_rho": float(r), "cgga_n": int(len(cg)),
               "cgga_events": int(cg.event.sum())}, f, ensure_ascii=False, indent=2)
print("\n完成 -> RI_axis_signature_genes.tsv / cgga_riaxis_scores.tsv")
