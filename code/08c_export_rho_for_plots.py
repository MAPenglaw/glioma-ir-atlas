#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""导出画图用数据：各事件类型 rho 分布 + SNRPC 散点（低内存版）"""
import gc
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output"
MIN_DETECT = 0.90


def barcode(c):
    return str(c).strip().strip('"')[:15]


def rank_rows(A):
    A = np.asarray(A, dtype=float)
    out = np.empty_like(A)
    for i in range(A.shape[0]):
        out[i] = rankdata(A[i])
    return out


def row_corr(A, B):
    a = np.nan_to_num(A - np.nanmean(A, axis=1, keepdims=True))
    b = np.nan_to_num(B - np.nanmean(B, axis=1, keepdims=True))
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))


print("[1] 读 PSI 元信息")
meta = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t", index_col=0,
                   compression="gzip", usecols=[0, 1, 4], low_memory=False)
print(f"    {len(meta):,} 事件")

print("[2] 读表达谱")
expr = pd.read_csv(f"{SRC}/GBM_HiSeqV2.gz", sep="\t", index_col=0, compression="gzip")
expr.columns = [barcode(c) for c in expr.columns]
expr = expr.loc[~expr.index.duplicated()]

print("[3] 逐块读 PSI 数值")
rows = []
reader = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t", index_col=0,
                     compression="gzip", chunksize=8000, low_memory=False)
keep = meta.index
for ci, ch in enumerate(reader):
    v = ch.iloc[:, 13:]
    v.columns = [barcode(c) for c in v.columns]
    tum = [c for c in v.columns if len(c) >= 15 and c[13:15] == "01"]
    v = v[tum].loc[:, ~v[tum].columns.duplicated()]
    v = v.apply(pd.to_numeric, errors="coerce")
    common = [c for c in v.columns if c in set(expr.columns)]
    v = v[common]
    rows.append(v.loc[v.index.intersection(keep)])
    del ch, v
    gc.collect()
P = pd.concat(rows)
del rows
gc.collect()
print(f"    PSI: {P.shape[0]:,} 事件 × {P.shape[1]} 样本")

detect = P.notna().sum(1).values / P.shape[1]
Psd = P.std(axis=1, skipna=True).values
gene = meta.loc[P.index, "GeneSymbol"].astype(str).values
st = meta.loc[P.index, "SpliceType"].astype(str).values
in_expr = np.isin(gene, np.asarray(expr.index))
base = (detect >= MIN_DETECT) & in_expr & (Psd > 0)
print(f"    可分析事件 {int(base.sum()):,}")

out = []
for t in ["IR", "ES", "A5"]:
    m = base & (st == t)
    Pk = P.loc[m, :]
    med = Pk.median(axis=1).values
    A = Pk.values.copy()
    nan_pos = np.isnan(A)
    idx = np.where(nan_pos)
    A[idx] = np.take(med, idx[0])
    Rp = rank_rows(A)
    Re = rank_rows(expr.loc[gene[m], Pk.columns].values)
    rho = row_corr(Rp, Re)
    rho = rho[~np.isnan(rho)]
    print(f"    {t}: n={len(rho):,} 负相关 {np.mean(rho < 0):.1%} 中位 {np.median(rho):+.4f}")
    out.append(pd.DataFrame({"type": t, "rho": rho}))
    del Pk, A, Rp, Re
    gc.collect()

pd.concat(out).to_csv(f"{OUT}/G1_rho_distribution.tsv", sep="\t", index=False)
print(f"[4] 已输出 {OUT}/G1_rho_distribution.tsv")

ann = pd.read_csv(f"{OUT}/sample_RI_group_annotated.tsv", sep="\t", index_col=0)
ann.index = [str(i)[:15] for i in ann.index]
cols = list(P.columns)
sc = pd.DataFrame({"sample": cols, "RI_burden": ann.loc[cols, "RI_burden"].values})
for g in ["SNRPC", "SNRPG", "SNRNP70", "ZMAT2", "ISY1"]:
    if g in expr.index:
        sc[g] = expr.loc[g, cols].values
sc.to_csv(f"{OUT}/G2_scatter_data.tsv", sep="\t", index=False)
print(f"[5] 已输出 {OUT}/G2_scatter_data.tsv")
