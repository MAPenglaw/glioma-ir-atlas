#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""30 — 全癌种 IR PSI 流式扫描（Pass 1）。

数据源：data/merge_graphs_intron_retention_C2.confirmed.txt.gz
        10,018 样本列 × 123,390 IR 事件（Kahles 2018 / SplAdder，GDC Open-Access）

目的：统计每个事件在「原发肿瘤（sampleTypeID=01）」中的检出样本数，
      为后续筛选「泛癌统一高检出事件集」提供依据。

为什么只统计 01：burden 只在原发肿瘤间可比；正常组织(11)/复发(02)/转移(06)
      样本深度与构成不同，混进来会污染检出率。

输出（output_pancan/）：
  col_meta.tsv            列索引 → sample_id / projectID / sampleTypeID
  event_detection.tsv     event_id + 在 01 样本中的检出数
  sample_detected.tsv     sample_id + 该样本检出的事件数
  is01.npy                01 列掩码
"""
import gzip
import csv
import os
import time
import numpy as np
import pandas as pd

DATA = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/data"
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output_pancan"
os.makedirs(OUT, exist_ok=True)

PSI = os.path.join(DATA, "merge_graphs_intron_retention_C2.confirmed.txt.gz")
SAMPLE_LIST = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
               "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/"
               "splicing_analysis/SplAdder_SampleList_TCGA.csv.gz")

N_META = 6          # event_id / event_type / event_chr / event_coordinates /
                    # alt_region_coordinates / gene_name
CHUNK = 2000        # 每块行数：2000 × 8300 float32 ≈ 66 MB

t0 = time.time()

# ---------------------------------------------------------------- 1. 样本映射
proj_of, type_of = {}, {}
with gzip.open(SAMPLE_LIST, "rt") as f:
    for row in csv.DictReader(f):
        sid = row["sampleID"].strip()
        proj_of[sid] = row["projectID"].strip()
        type_of[sid] = row["sampleTypeID"].strip()
print(f"[map] 样本清单 {len(proj_of)} 条", flush=True)

# ---------------------------------------------------------------- 2. 表头
with gzip.open(PSI, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
print(f"[hdr] 总列数 {len(header)} | 前6列 {header[:6]}", flush=True)

samp_cols = header[N_META:]
sids = ["-".join(c.split("-")[:4]) for c in samp_cols]
projs = [proj_of.get(s, "") for s in sids]
types = [type_of.get(s, "") for s in sids]
is01 = np.array([t == "01" for t in types])
print(f"[hdr] 样本列 {len(sids)} | 原发01 {int(is01.sum())} | "
      f"未匹配清单 {int(sum(1 for p in projs if p == ''))}", flush=True)

np.save(os.path.join(OUT, "is01.npy"), is01)
with open(os.path.join(OUT, "col_meta.tsv"), "w") as fh:
    fh.write("col_index\tsample_id\tprojectID\tsampleTypeID\n")
    for i, (s, p, t) in enumerate(zip(sids, projs, types)):
        fh.write(f"{i}\t{s}\t{p}\t{t}\n")

# ---------------------------------------------------------------- 3. 流式扫描
idx01 = np.where(is01)[0]
use_cols = [samp_cols[i] for i in idx01]
dtype_map = {c: "float32" for c in use_cols}
dtype_map[header[0]] = "str"

ev_det = []            # 每事件在 01 中的检出数
ev_ids = []            # event_id
samp_det = np.zeros(len(idx01), dtype=np.int64)   # 每样本检出事件数
n_rows = 0

reader = pd.read_csv(
    PSI, sep="\t", chunksize=CHUNK, usecols=[header[0]] + use_cols,
    dtype=dtype_map, low_memory=False, na_values=["", "NA", "nan", "NaN"])

for ci, chunk in enumerate(reader):
    ev_ids.extend(chunk[header[0]].tolist())
    vals = chunk[use_cols].to_numpy(dtype=np.float32, copy=False)
    det = np.count_nonzero(~np.isnan(vals), axis=1)
    ev_det.append(det)
    samp_det += np.count_nonzero(~np.isnan(vals), axis=0)
    n_rows += vals.shape[0]
    if (ci + 1) % 10 == 0:
        el = time.time() - t0
        print(f"[scan] chunk {ci+1} | 行 {n_rows} | "
              f"用时 {el/60:.1f} min", flush=True)

ev_det = np.concatenate(ev_det)
print(f"[scan] 完成：{n_rows} 事件 × {len(idx01)} 原发样本 | "
      f"总用时 {(time.time()-t0)/60:.1f} min", flush=True)

n01 = len(idx01)
with open(os.path.join(OUT, "event_detection.tsv"), "w") as fh:
    fh.write("event_id\tn_detected\tdetection_rate\n")
    for e, d in zip(ev_ids, ev_det):
        fh.write(f"{e}\t{d}\t{d/n01:.6f}\n")

sids01 = [sids[i] for i in idx01]
projs01 = [projs[i] for i in idx01]
with open(os.path.join(OUT, "sample_detected.tsv"), "w") as fh:
    fh.write("sample_id\tprojectID\tn_detected\n")
    for s, p, d in zip(sids01, projs01, samp_det):
        fh.write(f"{s}\t{p}\t{d}\n")

# ---------------------------------------------------------------- 4. 阈值分布
print("\n=== 泛癌(01)事件检出率分位数 ===", flush=True)
for th in [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99]:
    k = int((ev_det >= th * n01).sum())
    print(f"  检出率 ≥ {th:.2f} : {k} 个事件", flush=True)
print(f"  全检出(=1.00)  : {int((ev_det == n01).sum())} 个事件", flush=True)
print(f"\n中位数检出数 {np.median(ev_det):.0f} / {n01}", flush=True)

# 各癌种样本检出事件数（质控）
sd = pd.DataFrame({"projectID": projs01, "n_detected": samp_det})
print("\n=== 各癌种样本检出事件数（中位数） ===", flush=True)
g = sd.groupby("projectID")["n_detected"].agg(["count", "median"]).sort_values("count", ascending=False)
print(g.to_string(), flush=True)
print("\n[done]", flush=True)
