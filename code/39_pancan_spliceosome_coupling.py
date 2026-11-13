#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""39 — 泛癌剪接体耦合定位：32 癌种上「剪接体亚复合体 vs IR burden」。

核心问题：B 线发现 GBM 中 Sm 核心与 IR burden 负相关、LGG 中反转为正相关，
这个「GBM 特异的反转」到底是胶质瘤特例，还是泛癌坐标系里的一个极端位置？

设计要点
--------
1. 表达值在**每个癌种内部**做 z-score（跨癌种基线差异巨大，不可直接比）。
2. 亚复合体用成员基因 z 分数的均值（规避 7/7 伪复制）。
3. 相关用 Spearman，报告 ρ + 95% CI。
4. 阳性对照：Proliferation / Tumour-content markers（预期约等于 0）。

输出:
  output_pancan/pancan_sf_coupling.tsv      每个癌种 × 每个面板的 ρ/CI/p
  output_pancan/pancan_sf_expr.tsv          样本级面板分数（供绘图）
"""
import os
import gzip
import json
import sys
import time
import numpy as np
import pandas as pd
import urllib.request
from scipy import stats as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output_pancan")
EXPR = os.path.join(ROOT, "data", "pancan_expr",
                    "tcga_RSEM_gene_tpm.gz")
MAPF = os.path.join(OUT, "gene_ensg_map.tsv")

GROUPS = {
    "Sm_core": ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"],
    "SR_proteins": ["SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6",
                    "SRSF7", "SRSF8", "SRSF9", "SRSF10", "SRSF11", "SRSF12"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "AURKA",
                      "PLK1", "RRM2", "TYMS", "MCM2", "E2F1", "FOXM1"],
    "Tumor_content": ["GFAP", "SOX2", "OLIG2", "NES", "EGFR", "CD44",
                      "VIM", "PTPRZ1", "BCAN", "SPARC"],
    "NineteenComplex": ["ISY1", "ZMAT2", "CDC5L", "BCAS2", "BUD31", "CRNKL1",
                        "AQR", "XAB2", "PPIE", "SYF2", "CWC15", "CTNNBL1"],
}
SINGLE = ["SNRPC"]            # 主角单独看
ALL_GENES = sorted({g for v in GROUPS.values() for g in v} | set(SINGLE))


# ----------------------------------------------------------- 1. ENSG 映射
def build_map():
    # 增量缓存：每查到一个立即落盘，进程被杀也不丢（前台 >120s 会被 SIGKILL）
    done = {}
    if os.path.exists(MAPF):
        for _, r in pd.read_csv(MAPF, sep="\t", header=None,
                                names=["symbol", "ensg"],
                                dtype=str).iterrows():
            if r["ensg"]:
                done[r["symbol"]] = r["ensg"]
    todo = [g for g in ALL_GENES if g not in done]
    if not todo:
        return done
    op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with open(MAPF, "a") as fh:
        for g in todo:
            url = (f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{g}"
                   f"?content-type=application/json")
            e1 = ""
            for attempt in range(3):
                try:
                    r = json.load(op.open(urllib.request.Request(url), timeout=60))
                    ids = [x["id"] for x in r
                           if x.get("id", "").startswith("ENSG")]
                    e1 = sorted(set(ids))[0] if ids else ""
                    break
                except Exception as e:
                    print(f"  [warn] {g} attempt {attempt+1}: {e}", flush=True)
                    time.sleep(3)
            done[g] = e1
            fh.write(f"{g}\t{e1}\n")
            fh.flush()
            print(f"  {g:12s} -> {e1}", flush=True)
            time.sleep(0.12)
    return done


# -------------------------------------------- 2. 流式解压抽取目标基因行
def extract_expr(wanted_ensg):
    """wanted_ensg: set of ENSG (无版本号)。返回 DataFrame: index=ENSG, cols=barcode"""
    want = {e.split(".")[0] for e in wanted_ensg if e}
    rows = {}
    hdr = None
    t0 = time.time()
    with gzip.open(EXPR, "rt") as f:
        for i, line in enumerate(f):
            if i == 0:
                hdr = line.rstrip("\n").split("\t")
                continue
            gid = line.split("\t", 1)[0]
            base = gid.split(".")[0]
            if base in want and base not in rows:
                rows[base] = line.rstrip("\n").split("\t")[1:]
            if i % 10000 == 0:
                print(f"    扫描 {i} 行, 已命中 {len(rows)}/{len(want)} "
                      f"({time.time()-t0:.0f}s)", flush=True)
            if len(rows) == len(want):
                print(f"    全部命中 @ {i} 行", flush=True)
                break
    print(f"  [expr] 命中 {len(rows)} / {len(want)} 个基因", flush=True)
    df = pd.DataFrame.from_dict(rows, orient="index", columns=hdr[1:])
    return df


def main():
    if not os.path.exists(EXPR):
        sys.exit(f"[error] 表达矩阵未就绪: {EXPR}")

    # ---- 1. burden ----
    bdf = pd.read_csv(os.path.join(OUT, "pancan_burden.tsv"), sep="\t")
    if "patient" not in bdf.columns:
        bdf["patient"] = bdf["sample_id"].str[:15]
    print(f"[1] burden: {len(bdf)} 患者, {bdf.projectID.nunique()} 癌种", flush=True)

    # ---- 2. 映射 + 抽取 ----
    print("[2] 构建 ENSG 映射 …", flush=True)
    smap = build_map()
    ensg_of = {g: smap.get(g, "") for g in ALL_GENES}
    print("[3] 流式解压抽取目标基因 …", flush=True)
    E = extract_expr([v for v in ensg_of.values() if v])
    E.columns = [c[:15] for c in E.columns]           # patient 15 位
    E = E.loc[:, ~E.columns.duplicated()]
    E = E.astype(float)
    # Xena 用 log2(TPM+0.001)，-9.9658 为地板；低于 -9 视为未检出
    E = E.mask(E < -9.0)
    print(f"  [expr] 基因 {E.shape[0]} × 患者 {E.shape[1]}", flush=True)

    # ENSG -> symbol 反查
    rev = {v: k for k, v in ensg_of.items() if v}
    E.index = [rev.get(i, i) for i in E.index]
    print("  [expr] 基因符号:", sorted(E.index.tolist()), flush=True)

    # ---- 3. 每个癌种内部 z-score + 亚复合体均值 ----
    rows, sample_rows = [], []
    for proj, sub in bdf.groupby("projectID"):
        pats = sub["patient"].to_numpy()
        com = [p for p in pats if p in E.columns]
        if len(com) < 40:
            print(f"  [skip] {proj}: 表达交集 {len(com)} < 40", flush=True)
            continue
        Es = E[com]
        Z = Es.sub(Es.mean(1), axis=0).div(Es.std(1).replace(0, np.nan), axis=0)
        bv = sub.set_index("patient").loc[com, "RI_burden"].to_numpy()

        panels = {}
        for name, genes in GROUPS.items():
            have = [g for g in genes if g in Z.index]
            if len(have) >= 3:
                panels[name] = Z.loc[have].mean(0).to_numpy()
        for g in SINGLE:
            if g in Z.index:
                panels[g] = Z.loc[g].to_numpy()

        for name, vals in panels.items():
            ok = ~np.isnan(vals) & ~np.isnan(bv)
            if ok.sum() < 40:
                continue
            r = st.spearmanr(vals[ok], bv[ok])
            n = int(ok.sum())
            # Fisher z 95% CI
            z = np.arctanh(np.clip(r.statistic, -0.999, 0.999))
            se = 1 / np.sqrt(n - 3)
            lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
            rows.append({"projectID": proj, "panel": name, "n": n,
                         "rho": round(float(r.statistic), 4),
                         "ci_lo": round(float(lo), 4),
                         "ci_hi": round(float(hi), 4),
                         "p": float(r.pvalue)})
            for pv, bvv in zip(vals, bv):
                sample_rows.append({"patient": None, "projectID": proj,
                                    "panel": name, "score": pv,
                                    "RI_burden": bvv})
        print(f"  {proj}: n={len(com)}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "pancan_sf_coupling.tsv"), sep="\t",
               index=False)

    # ---- 4. 汇总 ----
    print("\n===== 泛癌剪接体耦合汇总（ρ with IR burden）=====", flush=True)
    for panel in ["Sm_core", "SR_proteins", "SNRPC", "Proliferation",
                  "Tumor_content", "NineteenComplex"]:
        s = res[res.panel == panel]
        if s.empty:
            continue
        print(f"\n{panel}  (n_cohorts={len(s)})", flush=True)
        print(f"  ρ 中位数 {s.rho.median():+.3f}  "
              f"范围 {s.rho.min():+.3f} ~ {s.rho.max():+.3f}", flush=True)
        print(f"  负相关癌种 {(s.rho < 0).sum()}/{len(s)}   "
              f"显著(q<0.05 前) p<0.05: {(s.p < 0.05).sum()}", flush=True)
        for _, r in s.sort_values("rho").head(4).iterrows():
            print(f"    最低 {r.projectID:12s} ρ={r.rho:+.3f} "
                  f"[{r.ci_lo:+.3f},{r.ci_hi:+.3f}] p={r.p:.2e}", flush=True)
        for _, r in s.sort_values("rho").tail(3).iterrows():
            print(f"    最高 {r.projectID:12s} ρ={r.rho:+.3f} "
                  f"[{r.ci_lo:+.3f},{r.ci_hi:+.3f}] p={r.p:.2e}", flush=True)
        gl = s[s.projectID.isin(["TCGA-GBM", "TCGA-LGG"])]
        for _, r in gl.iterrows():
            print(f"    ★ {r.projectID:12s} ρ={r.rho:+.3f} "
                  f"[{r.ci_lo:+.3f},{r.ci_hi:+.3f}] p={r.p:.2e}", flush=True)
    print("\n[done]")


if __name__ == "__main__":
    main()
