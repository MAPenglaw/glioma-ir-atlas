#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""31 — 提取「泛癌统一高检出事件」PSI 矩阵（Pass 2）。

在 30 号扫描的基础上，取在 8,301 例原发肿瘤中检出率 ≥ 90% 的 IR 事件
（与 B 线 GBM/LGG 的 burden 口径一致：检出 ≥90% 事件上的平均 PSI），
一次性抽出矩阵存盘，后续所有泛癌分析都基于它，不再重复读 937 MB 源文件。

输出（output_pancan/）：
  mat_hq.npy          (n_event × n_sample) float32，行为事件、列为原发肿瘤
  mat_hq_events.tsv   event_id / gene_name / 泛癌检出率
  mat_hq_samples.tsv  sample_id / projectID
"""
import os
import time
import numpy as np
import pandas as pd

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
DATA = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/data"
OUT = os.path.join(ROOT, "output_pancan")

PSI = os.path.join(DATA, "merge_graphs_intron_retention_C2.confirmed.txt.gz")
DET = os.path.join(OUT, "event_detection.tsv")
COLS = os.path.join(OUT, "col_meta.tsv")

N_META = 6
CHUNK = 3000
THRESH = 0.90

t0 = time.time()

# ------------------------------------------------- 1. 事件掩码（按文件行序）
det = pd.read_csv(DET, sep="\t")
mask = (det["detection_rate"].values >= THRESH)
n_keep = int(mask.sum())
print(f"[evt] 总事件 {len(det)} | 泛癌检出≥{THRESH:.2f} 保留 {n_keep}", flush=True)

row_of = np.full(len(mask), -1, dtype=np.int64)
row_of[mask] = np.arange(n_keep)

# ------------------------------------------------- 2. 01 样本列
# 注意：col_meta.tsv 的 sampleTypeID 是 "01" 这样的零填充字符串，
# 用默认 dtype 读会被 pandas 解析成整数 1，导致 == "01" 筛选为 0 条。
# 因此掩码直接复用 30 号保存的 is01.npy，名称列一律按字符串读。
is01 = np.load(os.path.join(OUT, "is01.npy"))
col = pd.read_csv(COLS, sep="\t", dtype=str, keep_default_na=False)
idx01 = np.where(is01)[0]
sids01 = col["sample_id"].values[is01]
projs01 = col["projectID"].values[is01]
n01 = len(idx01)
print(f"[col] 原发01 样本 {n01}", flush=True)

mat = np.empty((n_keep, n01), dtype=np.float32)
filled = np.zeros(n_keep, dtype=bool)

# ------------------------------------------------- 3. 表头与列名
with __import__("gzip").open(PSI, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
samp_cols = header[N_META:]
use_cols = [samp_cols[i] for i in idx01]
dtype_map = {c: "float32" for c in use_cols}
dtype_map[header[0]] = "str"

# ------------------------------------------------- 4. 流式填充
pos = 0
reader = pd.read_csv(
    PSI, sep="\t", chunksize=CHUNK, usecols=[header[0]] + use_cols,
    dtype=dtype_map, low_memory=False, na_values=["", "NA", "nan", "NaN"])

for ci, chunk in enumerate(reader):
    vals = chunk[use_cols].to_numpy(dtype=np.float32, copy=False)
    n = vals.shape[0]
    m = mask[pos:pos + n]
    if m.any():
        tgt = row_of[pos:pos + n][m]
        mat[tgt] = vals[m]
        filled[tgt] = True
    pos += n
    if (ci + 1) % 20 == 0:
        print(f"[ext] chunk {ci+1} | 已读行 {pos} | "
              f"用时 {(time.time()-t0)/60:.1f} min", flush=True)

assert filled.all(), f"有 {int((~filled).sum())} 个保留事件未被填充"
assert pos == len(det), f"行数不符：{pos} vs {len(det)}"
print(f"[ext] 矩阵填充完成 {mat.shape} | "
      f"用时 {(time.time()-t0)/60:.1f} min", flush=True)

# ------------------------------------------------- 5. 存盘
np.save(os.path.join(OUT, "mat_hq.npy"), mat)
print(f"[sav] mat_hq.npy {mat.nbytes/1e9:.2f} GB", flush=True)

ev = det.loc[mask, ["event_id", "detection_rate"]].copy()
ev["gene_name"] = np.nan
# gene_name 需要再读一次表头列（只取保留行）
ev_ids_keep = set(ev["event_id"].tolist())
gene_map = {}
reader2 = pd.read_csv(
    PSI, sep="\t", chunksize=CHUNK, usecols=[header[0], header[5]],
    dtype={header[0]: "str", header[5]: "str"}, low_memory=False)
for chunk in reader2:
    sub = chunk[chunk[header[0]].isin(ev_ids_keep)]
    for a, b in zip(sub[header[0]], sub[header[5]]):
        gene_map[a] = b
ev["gene_name"] = ev["event_id"].map(gene_map)
ev.to_csv(os.path.join(OUT, "mat_hq_events.tsv"), sep="\t", index=False)

pd.DataFrame({"sample_id": sids01, "projectID": projs01}).to_csv(
    os.path.join(OUT, "mat_hq_samples.tsv"), sep="\t", index=False)

print(f"[sav] mat_hq_events.tsv / mat_hq_samples.tsv", flush=True)
print(f"[nan] 矩阵缺失率 {np.isnan(mat).mean():.4f}", flush=True)
print("[done]", flush=True)
