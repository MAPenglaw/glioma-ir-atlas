#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 6 — GSE182109 单细胞定位：RI 轴签名基因在哪些细胞类型表达最高
细胞类型注释复用既有逻辑（PTPRC/CD3E/CD8A/AIF1/CD68/C1QA/SOX2/EGFR/CLDN5/COL1A1/PLP1）
"""
import gzip
import json
import os

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

SC = "/Users/peng/课题/省卫健委/sc_data/GSE182109"
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"

sig = pd.read_csv(f"{OUT}/output/RI_axis_signature_genes.tsv", sep="\t")
POS = sig.head(100)["gene"].tolist()
NEG = sig.tail(100)["gene"].tolist()
MARKERS = ["PTPRC", "CD3E", "CD8A", "AIF1", "CD68", "C1QA", "SOX2", "EGFR",
           "CLDN5", "COL1A1", "PLP1", "SNRPC", "UQCRC1"]
WANT = sorted(set(POS + NEG + MARKERS))


def load_sample(prefix, sample_id):
    feat = pd.read_csv(prefix + "_features.tsv.gz", sep="\t", header=None)
    syms = feat[1].astype(str).values
    pos = {}
    for i, s in enumerate(syms):
        if s in WANT and s not in pos:
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
    df.index = [b + "|" + sample_id for b in bc]
    df["total_counts"] = totals
    df["sample"] = sample_id
    return df


def main():
    prefixes = sorted(set(f[: -len("_barcodes.tsv.gz")] for f in os.listdir(SC)
                          if f.endswith("_barcodes.tsv.gz")))
    print(f"样本数: {len(prefixes)}", flush=True)
    dfs = []
    for n, p in enumerate(prefixes):
        sid = os.path.basename(p).split("_", 1)[1]
        try:
            d = load_sample(os.path.join(SC, p), sid)
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
    norm = np.log1p(ad[genes].div(ad.total_counts, axis=0) * 1e4)
    norm["sample"] = ad["sample"].values

    ptprc = norm["PTPRC"].values > 0
    tmark = (norm["CD3E"].values > 0) | (norm["CD8A"].values > 0)
    mmark = (norm["AIF1"].values > 0) | (norm["CD68"].values > 0) | (norm["C1QA"].values > 0)
    tum = (norm["SOX2"].values > 0) | (norm["EGFR"].values > 0)
    norm["cell_type"] = np.where(
        ptprc & tmark, "T_cell",
        np.where(ptprc & mmark, "Myeloid",
        np.where(ptprc, "Immune_other",
        np.where(tum, "Tumor",
        np.where(norm["CLDN5"].values > 0, "Endothelial",
        np.where(norm["COL1A1"].values > 0, "Fibroblast",
        np.where(norm["PLP1"].values > 0, "Oligo_like", "Other")))))))
    print("细胞类型分布:", norm["cell_type"].value_counts().to_dict(), flush=True)

    def mscore(cols):
        cs = [c for c in cols if c in norm.columns]
        return norm[cs].mean(axis=1).values if cs else np.full(len(norm), np.nan)

    norm["RIaxis_pos"] = mscore(POS)
    norm["RIaxis_neg"] = mscore(NEG)
    norm["RIaxis_score"] = norm["RIaxis_pos"] - norm["RIaxis_neg"]

    res = {}
    ct = norm.groupby("cell_type")[["RIaxis_pos", "RIaxis_neg", "RIaxis_score",
                                    "SNRPC", "UQCRC1"]].mean()
    cnt = norm.groupby("cell_type").size()
    ct["n_cells"] = cnt
    print("\n各细胞类型（平均 log1p-CP10K）:")
    print(ct.round(4).to_string())
    ct.to_csv(f"{OUT}/output/sc_celltype_RIaxis.tsv", sep="\t")

    res["n_cells"] = int(len(norm))
    res["cell_type_counts"] = norm["cell_type"].value_counts().to_dict()
    res["celltype_mean_score"] = ct["RIaxis_score"].round(5).to_dict()
    res["n_pos_genes"] = len([c for c in POS if c in norm.columns])
    res["n_neg_genes"] = len([c for c in NEG if c in norm.columns])

    # 肿瘤细胞 vs 非肿瘤的细胞级检验（近似，用均值差 + Mann-Whitney）
    from scipy.stats import mannwhitneyu
    a = norm.loc[norm.cell_type == "Tumor", "RIaxis_pos"]
    b = norm.loc[norm.cell_type != "Tumor", "RIaxis_pos"]
    if len(a) > 100 and len(b) > 100:
        u, p = mannwhitneyu(a, b, alternative="greater")
        res["tumor_vs_nontumor_pos_p"] = float(p)
        res["tumor_mean"] = float(a.mean())
        res["nontumor_mean"] = float(b.mean())
        print(f"\n恶性细胞 vs 非恶性: RI轴正相关基因均值 {a.mean():.4f} vs {b.mean():.4f}, "
              f"Mann-Whitney p={p:.2e}")

    norm[["cell_type", "sample", "RIaxis_pos", "RIaxis_neg", "RIaxis_score",
          "SNRPC"]].to_csv(f"{OUT}/output/sc_cell_RIaxis.tsv.gz",
                           sep="\t", compression="gzip", index=False)
    with open(f"{OUT}/output/step6_sc_summary.json", "w") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("\n完成 -> sc_celltype_RIaxis.tsv / step6_sc_summary.json")


if __name__ == "__main__":
    main()
