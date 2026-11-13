#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""17 — Kahles 2018 数据的机制验证：增殖/纯度解耦 + 剪接体核心关联。

用法: python 17_km_confound.py GBM | LGG
输出 ri_atlas/output_km/{cancer}/confound.tsv + sf_corr.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/07_原始项目数据_2026_splicing_analysis/"
       "2026 splicing analysis/splicing_analysis")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")

# 旧 HGNC 符号映射（Xena hg19 表达谱用 SFRS 而非 SRSF）
ALIAS = {"SRSF1": "SFRS1", "SRSF2": "SFRS2", "SRSF3": "SFRS3", "SRSF4": "SFRS4",
         "SRSF5": "SFRS5", "SRSF6": "SFRS6", "SRSF7": "SFRS7", "SRSF8": "SFRS8",
         "SRSF9": "SFRS9", "SRSF10": "SFRS13A", "SRSF11": "SFRS11", "SRSF12": "SFRS12"}

GROUPS = {
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "BUB1",
                      "BUB1B", "AURKA", "AURKB", "PLK1", "CENPA", "CENPF", "RRM2",
                      "TYMS", "MCM2", "MCM6", "E2F1", "FOXM1", "TTK"],
    "TumorPurity_marker": ["GFAP", "SOX2", "OLIG2", "NES", "EGFR", "CD44",
                           "VIM", "PTPRZ1", "BCAN", "SPARC"],
    "Neuronal_normal": ["RBFOX3", "SYT1", "SNAP25", "GRIN1", "GABRA1",
                        "MAP2", "NEFL", "NEFM", "STMN2", "SYP"],
    "Sm_core": ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"],
    "SR_proteins": ["SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6",
                    "SRSF7", "SRSF8", "SRSF9", "SRSF10", "SRSF11", "SRSF12"],
    "NineTeen_Bact": ["ISY1", "ZMAT2", "CDC5L", "BCAS2", "BUD31", "CRNKL1",
                      "AQR", "XAB2", "PPIE", "SYF2", "CWC15", "CTNNBL1"],
}


def barcode(c):
    return str(c).strip().strip('"')[:15]


def main():
    cancer = sys.argv[1] if len(sys.argv) > 1 else "GBM"
    outdir = os.path.join(OUT, cancer)
    burden = pd.read_csv(os.path.join(outdir, "burden.tsv"), sep="\t", index_col=0)["RI_burden"]
    burden.index = [barcode(i) for i in burden.index]
    # 口径对齐到表达谱（前 15 位），重复样本（01A/01B）取平均
    if burden.index.duplicated().any():
        burden = burden.groupby(level=0).mean()

    # 表达谱
    if cancer == "GBM":
        expr_path = os.path.join(SRC, "GBM_HiSeqV2.gz")
        comp = "gzip"
    else:
        expr_path = os.path.join(ROOT, "data", "LGG_HiSeqV2_full.txt")
        comp = None
    expr = pd.read_csv(expr_path, sep="\t", index_col=0, compression=comp, low_memory=False)
    expr.columns = [barcode(c) for c in expr.columns]
    expr = expr.loc[:, ~expr.columns.duplicated()]
    # 基因名去重（同符号取均值）
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    common = [c for c in burden.index if c in expr.columns]
    print(f"[1] {cancer}: burden 样本 {len(burden)}，表达谱交集 {len(common)}")
    E = expr[common]
    Z = E.sub(E.mean(1), axis=0).div(E.std(1).replace(0, np.nan), axis=0)

    # 模块分数
    scores = pd.DataFrame(index=common)
    for name, genes in GROUPS.items():
        have = [g for g in genes if g in Z.index]
        # 旧符号回退
        if len(have) < 3:
            have = [g for g in genes if ALIAS.get(g, g) in Z.index]
        if len(have) >= 3:
            idx = [g if g in Z.index else ALIAS.get(g, g) for g in have]
            idx = [g for g in idx if g in Z.index]
            scores[name] = Z.loc[idx].mean(axis=0)
        print(f"    {name:<18} {len(have)}/{len(genes)} 基因")

    # 单基因
    def single(g):
        if g in Z.index:
            return Z.loc[g]
        a = ALIAS.get(g, g)
        if a in Z.index:
            return Z.loc[a]
        return None

    print("\n[2] 混淆 + 机制 与 RI burden 的 Spearman 相关:")
    rows = []
    b = burden.reindex(common)
    for name in list(GROUPS.keys()):
        if name not in scores:
            continue
        r, p = st.spearmanr(b, scores[name])
        rows.append({"gene/panel": name, "rho": r, "p": p, "n": len(common)})
        print(f"    {name:<18} rho={r:+.3f}  p={p:.2e}")

    # 单基因：SNRPC 和 Sm 核心、SR 蛋白逐个
    single_genes = (["SNRPC"] + GROUPS["Sm_core"] +
                    ["SRSF1", "SRSF2", "SRSF3", "SRSF6", "SRSF7", "SRSF8"])
    print("\n[3] 单基因与 RI burden:")
    srows = []
    for g in single_genes:
        v = single(g)
        if v is None:
            continue
        v = v.reindex(common)
        r, p = st.spearmanr(b, v)
        srows.append({"gene": g, "rho": r, "p": p})
        print(f"    {g:<12} rho={r:+.3f}  p={p:.2e}")

    # 汇总 Sm 核心 / SR 蛋白方向
    sm = [r["rho"] for r in srows if r["gene"] in GROUPS["Sm_core"]]
    srs = [r["rho"] for r in srows if r["gene"].startswith("SRSF")]
    if sm:
        print(f"\n[4] Sm 核心 {sum(1 for x in sm if x<0)}/{len(sm)} 负（中位 {np.median(sm):+.3f}）")
    if srs:
        print(f"    SR 蛋白 {sum(1 for x in srs if x>0)}/{len(srs)} 正（中位 {np.median(srs):+.3f}）")

    pd.DataFrame(rows).to_csv(os.path.join(outdir, "confound.tsv"), sep="\t", index=False)
    pd.DataFrame(srows).to_csv(os.path.join(outdir, "sf_corr.tsv"), sep="\t", index=False)
    print("\n完成 ->", outdir)


if __name__ == "__main__":
    main()
