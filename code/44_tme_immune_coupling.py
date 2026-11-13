#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""44 — 肿瘤微环境（TME）耦合探查：免疫/基质面板 vs IR burden（GBM + LGG）。

背景：B 线手稿只排除了 proliferation / tumour-marker content / neuronal，
      **从未系统检验免疫与基质浸润**。而 Limitations 第六条却写了
      "argues that the bulk signal is not driven by microenvironment composition"
      ——这是用单细胞间接论证的，缺 bulk 层面的直接证据。本脚本补上这个洞。

设计（与 39 号保持一致，便于并列成图）
------------------------------------
1. 表达 = Toil 统一 RSEM（tcga_RSEM_gene_tpm.gz），log2(TPM+0.001)。
2. 在**癌种内部**做 z-score（GBM 内、LGG 内分别）。
3. 面板分数 = 成员基因 z 分数均值（规避伪复制）。
4. Spearman ρ + Fisher 95% CI + MDE(80% power)。
5. 必带对照：Sm_core（预期 GBM 负、LGG 正）、Proliferation（预期≈0）、
   Tumour_content（预期≈0）。对照跑不对说明管道有问题。

输出
----
  output_pancan/tme_immune_coupling.tsv   癌种 × 面板 的 ρ/CI/p/n_genes
  output_pancan/tme_sample_scores.tsv     样本级面板分数（供绘图）
"""
import os
import gzip
import json
import time
import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
from scipy import stats as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output_pancan")
EXPR = os.path.join(ROOT, "data", "pancan_expr", "tcga_RSEM_gene_tpm.gz")
MAPF = os.path.join(OUT, "gene_ensg_map.tsv")
BURDEN = {
    "GBM": os.path.join(ROOT, "output_km", "GBM", "burden.tsv"),
    "LGG": os.path.join(ROOT, "output_km", "LGG", "burden.tsv"),
}

TME_GROUPS = {
    "T_cell_general": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "CD8_T": ["CD8A", "CD8B", "GZMK", "CCL5"],
    "Treg": ["FOXP3", "IL2RA", "TIGIT", "CTLA4"],
    "NK": ["KLRD1", "NKG7", "GNLY", "KLRB1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B"],
    "Myeloid_general": ["LYZ", "ITGAM", "CD14", "FCGR3A", "AIF1", "CSF1R"],
    "Microglia": ["P2RY12", "CX3CR1", "TMEM119", "SALL1"],
    "M2_like": ["CD163", "MSR1", "MRC1", "CLEC7A"],
    "Neutrophil": ["FCGR3B", "S100A8", "S100A9", "CSF3R"],
    "DC": ["CD1C", "FCER1A", "CLEC10A"],
    "CAF_Stromal": ["FAP", "POSTN", "COL1A1", "COL1A2", "COL3A1",
                    "FN1", "DCN", "LUM", "ACTA2", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "TEK"],
    "IFN_gamma_Ayers": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1",
                        "IFNG", "ISG20", "GBP1", "CD274", "CXCR6"],
    "Cytotoxic": ["GZMA", "PRF1", "GZMB", "GZMH"],
    "MHC_I": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "NLRC5"],
    "MHC_II": ["HLA-DRA", "HLA-DPA1", "HLA-DMA", "CIITA"],
    "Checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
}
CONTROL_GROUPS = {
    "Sm_core": ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "AURKA",
                      "PLK1", "RRM2", "TYMS", "MCM2", "E2F1", "FOXM1"],
    "Tumor_content": ["GFAP", "SOX2", "OLIG2", "NES", "EGFR", "CD44",
                      "VIM", "PTPRZ1", "BCAN", "SPARC"],
}
GROUPS = {**TME_GROUPS, **CONTROL_GROUPS}
ALL_GENES = sorted({g for v in GROUPS.values() for g in v})


def build_map():
    """批量 POST 一次性拿全部基因；每个 symbol 保留全部候选 ENSG（逗号分隔），
    由 load_expr 后按平均表达量选主 ENSG（规避 HLA 多等位基因错配）。"""
    done = {}
    if os.path.exists(MAPF):
        for _, r in pd.read_csv(MAPF, sep="\t", header=None,
                                names=["symbol", "ensg"], dtype=str).iterrows():
            if isinstance(r["ensg"], str) and r["ensg"]:
                done[r["symbol"]] = r["ensg"]
    todo = [g for g in ALL_GENES if not done.get(g)]
    if not todo:
        return done
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    def q(g):
        u = ("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"
             f"{urllib.parse.quote(g)}?content-type=application/json")
        for _ in range(3):
            try:
                with opener.open(u, timeout=25) as r:
                    v = json.loads(r.read().decode())
                cand = sorted({x.get("id") for x in v
                               if str(x.get("id", "")).startswith("ENSG")})
                return g, ",".join(c for c in cand if c)
            except Exception:
                time.sleep(1.5)
        return g, ""

    from concurrent.futures import ThreadPoolExecutor
    got = {}
    with ThreadPoolExecutor(12) as ex:
        for i, (g, v) in enumerate(ex.map(q, todo), 1):
            got[g] = v
            if v:
                with open(MAPF, "a") as fh:
                    fh.write(f"{g}\t{v}\n")
            if i % 10 == 0:
                print(f"  [map] {i}/{len(todo)}", flush=True)
    done.update(got)
    return done


def load_expr(ensg_set):
    """流式扫描 740MB：只对目标基因行做完整 split。"""
    rows = {}
    cols = None
    t0 = time.time()
    with gzip.open(EXPR, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = np.array(header[1:])
        for line in f:
            gid = line.split("\t", 1)[0].split(".")[0]
            if gid in ensg_set:
                vals = np.fromstring(line.rstrip("\n").split("\t", 1)[1],
                                     sep="\t", dtype=np.float32)
                rows[gid] = vals
    print(f"  [expr] {len(rows)}/{len(ensg_set)} genes, "
          f"{len(cols)} samples, {time.time()-t0:.0f}s")
    return pd.DataFrame(rows, index=cols)


def fisher_ci(r, n):
    if n < 10 or not np.isfinite(r) or abs(r) >= 0.999:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    return (np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se))


def mde(n, alpha=0.05, power=0.8):
    """80% 效能下可检出的最小 |rho|（近似）。"""
    za = st.norm.ppf(1 - alpha / 2)
    zb = st.norm.ppf(power)
    return float(np.tanh((za + zb) / np.sqrt(max(n - 3, 1))))


def main():
    os.makedirs(OUT, exist_ok=True)
    gmap = build_map()
    miss = [g for g in ALL_GENES if not gmap.get(g)]
    if miss:
        print(f"  [warn] unmapped (will be skipped): {miss}")
    cands = {g: [x for x in gmap.get(g, "").split(",") if x]
             for g in ALL_GENES}
    ensg_set = {x for v in cands.values() for x in v}
    expr = load_expr(ensg_set)
    # 每个 symbol 只保留平均表达最高的 ENSG（HLA 等多等位基因家族的主基因）
    chosen = {}
    for g, cs in cands.items():
        ok = [c for c in cs if c in expr.columns]
        if not ok:
            chosen[g] = None
            continue
        chosen[g] = max(ok, key=lambda c: float(expr[c].mean()))
    drop = [c for c in expr.columns
            if c not in set(v for v in chosen.values() if v)]
    expr = expr.drop(columns=drop)
    expr.columns = [next(g for g, v in chosen.items() if v == c)
                    for c in expr.columns]
    print(f"  [map] resolved {sum(1 for v in chosen.values() if v)}"
          f"/{len(ALL_GENES)} symbols")

    # 只保留原发肿瘤(-01)，再折叠到 patient key；否则同一 patient 的
    # -01/-11 会同时命中，造成行数多于 burden（此前 156 vs 151 的广播错误）
    keep = [i for i, s in enumerate(expr.index) if s.endswith("-01")]
    expr = expr.iloc[keep]
    expr.index = [s[:12] for s in expr.index]
    expr = expr[~expr.index.duplicated(keep="first")]

    out_rows = []
    sample_rows = []
    for cohort, bpath in BURDEN.items():
        bd = pd.read_csv(bpath, sep="\t", index_col=0)
        bd.index = [s[:12] for s in bd.index]
        bd = bd[~bd.index.duplicated(keep="first")]
        common = bd.index.intersection(expr.index)
        print(f"  [{cohort}] burden n={len(bd)}, matched n={len(common)}")
        e = expr.loc[common]
        y = bd.loc[common, "RI_burden"].values
        for gname, genes in GROUPS.items():
            es = [g for g in genes if g in e.columns]
            if len(es) < 2:
                continue
            sub = e[es].astype(float)
            z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
            score = z.mean(axis=1).values
            ok = np.isfinite(score) & np.isfinite(y)
            r, p = st.spearmanr(score[ok], y[ok])
            lo, hi = fisher_ci(r, int(ok.sum()))
            out_rows.append(dict(cohort=cohort, panel=gname,
                                 rho=r, ci_lo=lo, ci_hi=hi, p=p,
                                 n=int(ok.sum()), n_genes=len(es),
                                 mde80=mde(int(ok.sum()))))
            sample_rows.append(pd.DataFrame(
                {"cohort": cohort, "panel": gname, "sample": common,
                 "score": score, "RI_burden": y}))
    res = pd.DataFrame(out_rows)
    res.to_csv(os.path.join(OUT, "tme_immune_coupling.tsv"),
               sep="\t", index=False)
    pd.concat(sample_rows).to_csv(os.path.join(OUT, "tme_sample_scores.tsv"),
                                  sep="\t", index=False)
    print("\n=== TME immune coupling ===")
    print(res.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))


if __name__ == "__main__":
    main()
