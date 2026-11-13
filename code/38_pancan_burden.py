#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""38 — 导出泛癌样本级 IR burden（口径与 B 线 16 号一致）。

burden = 该癌种内检出率 ≥90% 的 IR 事件上的平均 PSI

输出:
  output_pancan/pancan_burden.tsv       样本级 (patient 15 位 barcode)
  output_pancan/pancan_burden_cohort.tsv 各癌种事件数与 burden 分布
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output_pancan")

M = np.load(os.path.join(OUT, "mat_hq.npy"), mmap_mode="r")
samples = pd.read_csv(os.path.join(OUT, "mat_hq_samples.tsv"), sep="\t")
ev = pd.read_csv(os.path.join(OUT, "mat_hq_events.tsv"), sep="\t")
print(f"[load] matrix {M.shape}, samples {len(samples)}", flush=True)

assert M.shape[1] == len(samples), "样本列不匹配"
assert M.shape[0] == len(ev), "事件行不匹配"

rows = []
for proj, sub in samples.groupby("projectID"):
    idx = sub.index.to_numpy()
    Ms = np.asarray(M[:, idx])
    det = (~np.isnan(Ms)).mean(axis=1)
    m_b = det >= 0.90
    n_b = int(m_b.sum())
    if n_b < 100:
        print(f"  [skip] {proj}: 高检出事件 {n_b} < 100", flush=True)
        continue
    b = np.nanmean(Ms[m_b, :], axis=0)
    # 高变异事件用于后续 PCA（此处只需 burden，但一并记录）
    sd = np.nanstd(Ms[m_b, :], axis=1)
    n_v = int((sd >= 0.10).sum())
    for sid, bv in zip(sub["sample_id"].to_numpy(), b):
        rows.append({"sample_id": sid, "projectID": proj,
                     "RI_burden": bv, "n_burden_events": n_b})
    print(f"  {proj}: n={len(idx)} events={n_b} (hivar {n_v}) "
          f"burden median={np.nanmedian(b):.4f}", flush=True)

df = pd.DataFrame(rows)
df["patient"] = df["sample_id"].str[:15]
# 同一 patient 多个样本（01A/01B）取平均
dup = df["patient"].duplicated().sum()
if dup:
    print(f"[dedup] {dup} 个重复 patient 样本，取均值", flush=True)
    df = (df.groupby(["patient", "projectID"], as_index=False)["RI_burden"]
            .mean())

df.to_csv(os.path.join(OUT, "pancan_burden.tsv"), sep="\t", index=False)
print(f"[save] pancan_burden.tsv  n={len(df)}", flush=True)

co = (df.groupby("projectID")["RI_burden"]
        .agg(n="size", median="median", mean="mean", sd="std")
        .reset_index().sort_values("n", ascending=False))
co.to_csv(os.path.join(OUT, "pancan_burden_cohort.tsv"), sep="\t", index=False)
print(co.to_string(index=False), flush=True)
print("[done]")
