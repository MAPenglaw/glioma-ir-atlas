#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""21 — 把 G1（PSI–表达负耦合）迁移到 Kahles 2018 事件集。

与旧版 08 脚本的差异：
  - IR/ES 事件来自 Kahles 2018（gene_name 列是 ENSG ID）
  - ENSG → 新 symbol（gtf_gene_info.json）→ 旧 symbol（表达谱用旧 HGNC，SRSF→SFRS）
  - 主指标：负相关事件比例（NMD 预期 >50%），IR vs ES 配对基因内对比
"""
import os, json, gzip, csv
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
DATA = os.path.join(ROOT, "ri_atlas", "data")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
os.makedirs(OUT, exist_ok=True)
SRC = "/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis"
EXPR = os.path.join(SRC, "GBM_HiSeqV2.gz")
GTF = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc-as-glioma/cgga325/A_layer_n40/gtf_gene_info.json"
IR_PSI = os.path.join(DATA, "KM_IR_GBM.psi.gz")
IR_RAW = os.path.join(ROOT, "data", "merge_graphs_intron_retention_C2.confirmed.txt.gz")
ES_RAW = os.path.join(ROOT, "data", "merge_graphs_exon_skip_C2.confirmed.txt.gz")
SAMPLE_LIST = os.path.join(SRC, "SplAdder_SampleList_TCGA.csv.gz")

# SRSF(新) -> SFRS(旧) 及其他常见改名
ALIAS = {"SRSF1": "SFRS1", "SRSF2": "SFRS2", "SRSF3": "SFRS3", "SRSF4": "SFRS4",
         "SRSF5": "SFRS5", "SRSF6": "SFRS6", "SRSF7": "SFRS7", "SRSF8": "SFRS8",
         "SRSF9": "SFRS9", "SRSF10": "SFRS13A", "SRSF11": "SFRS11", "SRSF12": "SFRS12"}


def barcode(c):
    return str(c).strip().strip('"')[:15]


def ensg_base(e):
    """ENSG00000187634.6 -> ENSG00000187634"""
    return str(e).split(".")[0]


# ---------------------------------------------------------------- 1. 表达谱 + 映射
print("[1] 读表达谱 + ENSG→symbol 映射 ...", flush=True)
expr = pd.read_csv(EXPR, sep="\t", index_col=0, low_memory=False)
if expr.index.duplicated().any():
    expr = expr.groupby(level=0).mean()
expr.columns = [barcode(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()]
expr_genes = set(expr.index)
print(f"    表达谱 {expr.shape[0]} 基因 × {expr.shape[1]} 样本")

gtf = json.load(open(GTF))
ensg2sym = {}
for e, v in gtf.items():
    if v and v[0]:
        ensg2sym[ensg_base(e)] = v[0]
print(f"    ENSG→symbol 映射 {len(ensg2sym)} 基因")

# 表达谱是旧 symbol，需要新 symbol -> 旧 symbol 的 alias；也直接尝试新 symbol
def resolve_sym(sym):
    if sym in expr_genes:
        return sym
    if sym in ALIAS and ALIAS[sym] in expr_genes:
        return ALIAS[sym]
    return None


# ---------------------------------------------------------------- 2. IR 事件 PSI + gene_name
print("[2] 读 Kahles IR PSI + gene_name ...", flush=True)
psi = pd.read_csv(IR_PSI, sep="\t", compression="gzip", low_memory=False)
psi_meta = psi[["event_id", "event_chr", "event_coordinates", "alt_region_coordinates"]]
psi_vals = psi.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
psi_vals.index = psi["event_id"].values
psi_vals.columns = [barcode(c) for c in psi_vals.columns]
psi_vals = psi_vals.loc[:, ~psi_vals.columns.duplicated()]

# gene_name 从原始文件流式读
gene_map = {}
with gzip.open(IR_RAW, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        gene_map[p[0]] = ensg_base(p[5])
print(f"    IR 事件 {len(psi_vals)}，gene_name 映射 {len(gene_map)}")


def compute_negfrac(vals, gene_map):
    """对每个事件，跨样本算 Spearman(PSI, 宿主基因表达)，返回负相关比例与中位 rho。"""
    # 样本对齐
    common = list(expr.columns.intersection(vals.columns))
    vals = vals[common]
    # 只保留有基因映射且基因在表达谱里的事件
    keep_rows = []
    keep_sym = []
    for i, eid in enumerate(vals.index):
        ensg = gene_map.get(eid, "")
        sym = ensg2sym.get(ensg, "")
        rsym = resolve_sym(sym) if sym else None
        if rsym:
            keep_rows.append(i)
            keep_sym.append(rsym)
    V = vals.iloc[keep_rows].values.astype(float)
    print(f"    可分析事件（有宿主基因）: {len(keep_rows):,}")
    # 宿主基因表达矩阵（事件 × 样本）
    E = expr.loc[keep_sym, common].values.astype(float)
    # Spearman = Pearson on ranks（逐事件）
    def rank_rows(M):
        return np.apply_along_axis(lambda r: stats.rankdata(r, nan_policy="omit"), 1, M)
    # 只对非全 NaN 的事件算
    rho = np.full(len(keep_rows), np.nan)
    for i in range(len(keep_rows)):
        p = V[i]; e = E[i]
        m = ~(np.isnan(p) | np.isnan(e))
        if m.sum() < 10:
            continue
        rho[i] = stats.spearmanr(p[m], e[m]).correlation
    neg = rho < 0
    frac = neg[~np.isnan(rho)].mean()
    med = np.nanmedian(rho)
    n_valid = (~np.isnan(rho)).sum()
    return frac, med, n_valid, rho


print("[3] 计算 IR 负相关比例 ...", flush=True)
ir_frac, ir_med, ir_n, ir_rho = compute_negfrac(psi_vals, gene_map)
print(f"    IR：负相关比例 {ir_frac:.3f}（{ir_n:,} 事件），中位 rho {ir_med:.4f}")
# 二项检验 vs 50%
p_binom = stats.binomtest(int(ir_frac * ir_n), ir_n, 0.5, alternative="greater").pvalue
print(f"    二项检验 vs 50%：p = {p_binom:.3g}")

# ---------------------------------------------------------------- 4. ES 事件（阴性对照）
print("[4] 读 Kahles ES 事件 ...", flush=True)
if os.path.exists(ES_RAW) and os.path.getsize(ES_RAW) > 500_000_000:
    # 流式读 ES 原始文件，只取 GBM 样本列 + gene_name
    with gzip.open(ES_RAW, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    meta_cols = header[:6]
    sample_cols = header[6:]
    # 用 sampleList 的 projectID + sampleTypeID 识别 GBM 原发肿瘤样本
    gbm_ids = set()
    with gzip.open(SAMPLE_LIST, "rt") as f:
        for row in csv.DictReader(f):
            if row.get("projectID", "").strip() == "TCGA-GBM" and row.get("sampleTypeID", "").strip() == "01":
                gbm_ids.add(row.get("sampleID", "").strip())
    def sample_id(col):
        return "-".join(col.split("-")[:4])
    gbm_idx = [i for i, c in enumerate(sample_cols) if sample_id(c) in gbm_ids]
    es_gene = {}
    es_rows = []
    with gzip.open(ES_RAW, "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            eid = p[0]
            ensg = ensg_base(p[5]) if len(p) > 5 else ""
            es_gene[eid] = ensg
            es_rows.append([p[0]] + [p[6 + i] for i in gbm_idx])
    es_vals = pd.DataFrame([r[1:] for r in es_rows], index=[r[0] for r in es_rows],
                           columns=[barcode(sample_cols[i]) for i in gbm_idx])
    es_vals = es_vals.loc[:, ~es_vals.columns.duplicated()]
    es_vals = es_vals.apply(pd.to_numeric, errors="coerce")
    print(f"    ES 事件 {len(es_vals)} × GBM 样本 {len(gbm_idx)}")
    es_frac, es_med, es_n, es_rho = compute_negfrac(es_vals, es_gene)
    print(f"    ES：负相关比例 {es_frac:.3f}（{es_n:,} 事件），中位 rho {es_med:.4f}")
    es_p_binom = stats.binomtest(int(es_frac * es_n), es_n, 0.5, alternative="greater").pvalue
    print(f"    ES 二项检验 vs 50%：p = {es_p_binom:.3g}")
else:
    print("    ES 文件未下载完成，跳过 ES 对照（ES_RAW 大小", os.path.getsize(ES_RAW) if os.path.exists(ES_RAW) else 0, "）")

print("\nDONE")
