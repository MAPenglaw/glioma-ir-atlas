#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""23 — 为图件重制导出绘图数据（Kahles 2018 双谱系口径）。

补齐 Fig1（逐事件检出率/变异度）、Fig3（PC1–PC10 方差 + PC 得分）、
Fig4（逐样本 burden + SNRPC/Sm核心/SR 蛋白表达 z）三张图所需数据。

用法: python 23_export_plot_data.py
输出到 ri_atlas/output_km/ ：
  fig1_event_stats.tsv   事件 × (cancer, detect, sd, meanPSI, hivar)
  fig3_pcevr.tsv         PC1–PC10 方差占比
  fig3_scores.tsv        sample, PC1, PC2, PC3, RI_burden, cancer
  fig4_expr.tsv          sample, RI_burden, SNRPC_z, Sm_core_z, SR_proteins_z, cancer
  fig1_summary.tsv       图谱汇总（事件数/检出数/高变异数/burden）
"""
import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/07_原始项目数据_2026_splicing_analysis/"
       "2026 splicing analysis/splicing_analysis")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
MIN_DETECT = 0.90
MIN_STD = 0.10

ALIAS = {"SRSF1": "SFRS1", "SRSF2": "SFRS2", "SRSF3": "SFRS3", "SRSF4": "SFRS4",
         "SRSF5": "SFRS5", "SRSF6": "SFRS6", "SRSF7": "SFRS7", "SRSF8": "SFRS8",
         "SRSF9": "SFRS9", "SRSF10": "SFRS13A", "SRSF11": "SFRS11", "SRSF12": "SFRS12"}
SM = ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"]
SR = ["SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6",
      "SRSF7", "SRSF8", "SRSF9", "SRSF10", "SRSF11", "SRSF12"]


def barcode(c):
    return str(c).strip().strip('"')[:15]


def zrow(x):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)


def main():
    ev_rows, evr_rows, sc_rows, ex_rows, sum_rows = [], [], [], [], []

    for cancer in ("GBM", "LGG"):
        psi_path = os.path.join(ROOT, "ri_atlas", "data", f"KM_IR_{cancer}.psi.gz")
        psi = pd.read_csv(psi_path, sep="\t", index_col=0, compression="gzip")
        vals = psi.iloc[:, 4:]
        vals = vals.apply(pd.to_numeric, errors="coerce")
        vals.columns = [str(c).split(".")[0] for c in vals.columns]
        vals = vals.loc[:, ~vals.columns.duplicated()]
        print(f"[{cancer}] 矩阵 {vals.shape[0]} × {vals.shape[1]}", flush=True)

        # ---- Fig1: 逐事件检出率 / 变异度 ----
        detect = vals.notna().sum(axis=1).values / vals.shape[1]
        sd = vals.std(axis=1, skipna=True).values
        mean_psi = vals.mean(axis=1, skipna=True).values
        ri_mask = detect >= MIN_DETECT
        hi = ri_mask & (sd >= MIN_STD)
        ev_rows.append(pd.DataFrame({
            "event": vals.index.astype(str), "cancer": cancer,
            "detect": detect, "sd": sd, "meanPSI": mean_psi,
            "retained": ri_mask, "hivar": hi}))
        print(f"  检出≥{MIN_DETECT}: {int(ri_mask.sum())}；高变异(std≥{MIN_STD}): {int(hi.sum())}",
              flush=True)

        # ---- burden ----
        burden = vals.loc[ri_mask, :].mean(axis=0, skipna=True)
        burden.name = "RI_burden"

        # ---- Fig3: PCA ----
        X = vals.loc[hi, :].T
        Xz = ((X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)).fillna(0.0)
        ncomp = min(10, Xz.shape[0] - 1)
        pca = PCA(n_components=ncomp)
        pc = pca.fit_transform(Xz.values)
        evr = pca.explained_variance_ratio_
        for i, v in enumerate(evr, start=1):
            evr_rows.append({"cancer": cancer, "PC": f"PC{i}", "var": float(v)})
        sc = pd.DataFrame({"sample": X.index.astype(str),
                           "PC1": pc[:, 0], "PC2": pc[:, 1], "PC3": pc[:, 2],
                           "RI_burden": burden.reindex(X.index).values,
                           "cancer": cancer})
        sc_rows.append(sc)
        print(f"  PC1 {evr[0]*100:.1f}% PC2 {evr[1]*100:.1f}%", flush=True)

        sum_rows.append({"cancer": cancer, "n_samples": int(vals.shape[1]),
                         "n_events": int(vals.shape[0]),
                         "n_detect": int(ri_mask.sum()), "n_hivar": int(hi.sum()),
                         "burden_mean": float(burden.mean()),
                         "burden_min": float(burden.min()),
                         "burden_max": float(burden.max())})

        # ---- Fig4: 逐样本表达 ----
        b = burden.copy()
        b.index = [barcode(i) for i in b.index]
        if b.index.duplicated().any():
            b = b.groupby(level=0).mean()
        if cancer == "GBM":
            expr = pd.read_csv(os.path.join(SRC, "GBM_HiSeqV2.gz"), sep="\t",
                               index_col=0, compression="gzip", low_memory=False)
        else:
            expr = pd.read_csv(os.path.join(ROOT, "data", "LGG_HiSeqV2_full.txt"),
                               sep="\t", index_col=0, low_memory=False)
        expr.columns = [barcode(c) for c in expr.columns]
        expr = expr.loc[:, ~expr.columns.duplicated()]
        if expr.index.duplicated().any():
            expr = expr.groupby(level=0).mean()
        common = [c for c in b.index if c in expr.columns]
        E = expr[common]
        Z = E.sub(E.mean(1), axis=0).div(E.std(1).replace(0, np.nan), axis=0)

        def score(genes):
            g = [ALIAS.get(x, x) for x in genes]
            g = [x for x in g if x in Z.index]
            return Z.loc[g].mean(axis=0, skipna=True)

        snrpc = Z.loc["SNRPC"] if "SNRPC" in Z.index else Z.loc[ALIAS.get("SNRPC", "SNRPC")]
        ex_rows.append(pd.DataFrame({
            "sample": common, "cancer": cancer,
            "RI_burden": b.reindex(common).values,
            "SNRPC_z": snrpc.reindex(common).values,
            "Sm_core_z": score(SM).reindex(common).values,
            "SR_proteins_z": score(SR).reindex(common).values}))
        print(f"  表达交集 n = {len(common)}", flush=True)

    pd.concat(ev_rows, ignore_index=True).to_csv(
        os.path.join(OUT, "fig1_event_stats.tsv"), sep="\t", index=False)
    pd.DataFrame(evr_rows).to_csv(os.path.join(OUT, "fig3_pcevr.tsv"), sep="\t", index=False)
    pd.concat(sc_rows, ignore_index=True).to_csv(
        os.path.join(OUT, "fig3_scores.tsv"), sep="\t", index=False)
    pd.concat(ex_rows, ignore_index=True).to_csv(
        os.path.join(OUT, "fig4_expr.tsv"), sep="\t", index=False)
    pd.DataFrame(sum_rows).to_csv(os.path.join(OUT, "fig1_summary.tsv"), sep="\t", index=False)
    print("DONE -> output_km/fig{1,3,4}_*.tsv")


if __name__ == "__main__":
    main()
