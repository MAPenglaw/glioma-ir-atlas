#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""15 — 从 Kahles 2018 全癌种 IR PSI 提取 GBM / LGG 子集。

数据源：data/merge_graphs_intron_retention_C2.confirmed.txt.gz
（GDC Open-Access，Cancer Cell 2018，SplAdder 全癌种 IR PSI）
样本映射：SplAdder_SampleList_TCGA.csv.gz（sampleID -> projectID / sampleTypeID）

输出（gzip TSV，5 列 meta + 原发肿瘤 01 样本列）：
  ri_atlas/data/KM_IR_GBM.psi.gz   (TCGA-GBM)
  ri_atlas/data/KM_IR_LGG.psi.gz   (TCGA-LGG)
"""
import gzip
import csv
import os

SRC = "/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis"
DATA = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/data"
OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/data"
os.makedirs(OUT, exist_ok=True)

PSI_FILE = os.path.join(DATA, "merge_graphs_intron_retention_C2.confirmed.txt.gz")
SAMPLE_LIST = os.path.join(SRC, "SplAdder_SampleList_TCGA.csv.gz")

N_META = 5  # event_id, event_type, event_chr, event_coordinates, alt_region_coordinates


def sample_id(col):
    """样本列名 TCGA-XX-XXXX-01A-...UUID.aligned -> barcode 前 16 位（TCGA-XX-XXXX-01A）"""
    return col.split("-")[0:4][0] if False else "-".join(col.split("-")[:4])


def barcode15(sid):
    """TCGA-XX-XXXX-01A -> 保留前 16 位（01A/01B 是不同样本，不能截成同一个 01）。"""
    return sid[:16]


def main():
    # 1. 读样本列表：sampleID -> (projectID, sampleTypeID)
    target = {"TCGA-GBM": set(), "TCGA-LGG": set()}
    with gzip.open(SAMPLE_LIST, "rt") as f:
        for row in csv.DictReader(f):
            sid = row.get("sampleID", "").strip()
            proj = row.get("projectID", "").strip()
            stype = row.get("sampleTypeID", "").strip()
            if proj in target and stype == "01":
                target[proj].add(sid)

    print("GBM 原发肿瘤:", len(target["TCGA-GBM"]), "| LGG 原发肿瘤:", len(target["TCGA-LGG"]))

    # 2. 读 PSI 表头，定位 GBM/LGG 样本列
    with gzip.open(PSI_FILE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    meta_cols = header[:N_META]
    samp_cols = header[N_META:]

    gbm_idx, lgg_idx = [], []
    gbm_names, lgg_names = [], []
    for i, col in enumerate(samp_cols):
        sid = sample_id(col)
        if sid in target["TCGA-GBM"]:
            gbm_idx.append(i)
            gbm_names.append(barcode15(sid))
        elif sid in target["TCGA-LGG"]:
            lgg_idx.append(i)
            lgg_names.append(barcode15(sid))

    print("命中 GBM 列:", len(gbm_idx), "| LGG 列:", len(lgg_idx))

    # 3. 流式写两个子集
    def write_subset(idx, names, outpath):
        w = gzip.open(outpath, "wt")
        w.write("\t".join(meta_cols + names) + "\n")
        n = 0
        with gzip.open(PSI_FILE, "rt") as f:
            f.readline()  # skip header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                meta = parts[:N_META]
                vals = [parts[N_META + i] for i in idx]
                w.write("\t".join(meta + vals) + "\n")
                n += 1
        w.close()
        print("写 %s: %d 事件" % (outpath, n))

    write_subset(gbm_idx, gbm_names, os.path.join(OUT, "KM_IR_GBM.psi.gz"))
    write_subset(lgg_idx, lgg_names, os.path.join(OUT, "KM_IR_LGG.psi.gz"))
    print("完成")


if __name__ == "__main__":
    main()
