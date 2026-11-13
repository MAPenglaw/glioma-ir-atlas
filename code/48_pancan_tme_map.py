#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""48 — 泛癌 IR burden × 微环境耦合地图（方案 C 第二部分，Fig9 的姊妹图）

问题
----
45 号在胶质瘤里发现：IR burden 与免疫/基质浸润负相关，且校正总表达量后**变强**
（suppressor 效应）。这引出两个泛癌层面的问题：

  Q1 这种「burden 越高、微环境越冷」的关系是胶质瘤特有，还是泛癌通则？
  Q2 GBM/LGG 在泛癌坐标里处在什么位置 —— 极端还是普通？

如果 Q1 的答案是「泛癌普遍」，那么胶质瘤的发现就不是一个孤立的队列巧合，
而是一个可推广的规律；GBM 在 Fig9 里已是 Sm core 耦合的泛癌极端，
若它在免疫轴上也处在极端，两图互为印证。

方法
----
- 31 个 TCGA 癌种（n >= 40），样本为原发肿瘤。
- 表达：Xena Toil tcga_RSEM_gene_tpm（log2(TPM+0.001)）。
  一次扫描同时取：面板基因行 + 每样本全局表达统计（TotalExpr）。
- 面板分数 = 面板内基因表达均值；癌种内做**秩偏相关**（k=0 时精确等于 Spearman）。
- 两个模型：
    M0  未校正
    M1  校正 TotalExpr（45 号已证明这是最强混杂，必须校正）
- 多重检验：每个面板内跨癌种 BH（31 个检验）。

输出
----
  output_pancan/pancan_tme_map.tsv        癌种 × 面板 的 rho/p/q/n（M0、M1）
  output_pancan/pancan_tme_samples.tsv    样本级（供绘图）
"""
import gzip
import os
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output_pancan")
EXPR = os.path.join(ROOT, "data", "pancan_expr", "tcga_RSEM_gene_tpm.gz")
MAPF = os.path.join(OUT, "gene_ensg_map.tsv")
BURDEN = os.path.join(OUT, "pancan_burden.tsv")
MIN_N = 40

# 泛癌地图聚焦的面板（17 个 TME 面板全做太密，选有代表性/可解释的）
PANELS = {
    "T_cell_general": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "CD8_T": ["CD8A", "CD8B", "GZMK", "CCL5"],
    "Myeloid_general": ["LYZ", "ITGAM", "CD14", "FCGR3A", "AIF1", "CSF1R"],
    "IFN_gamma_Ayers": ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1",
                        "IFNG", "ISG20", "GBP1", "CD274", "CXCR6"],
    "Cytotoxic": ["GZMA", "PRF1", "GZMB", "GZMH"],
    "CAF_Stromal": ["FAP", "POSTN", "COL1A1", "COL1A2", "COL3A1",
                    "FN1", "DCN", "LUM", "ACTA2", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "TEK"],
    # 主线对照：剪接体 Sm 核心（Fig9 同一对照，保证两图可比）
    "Sm_core": ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG"],
}
CONTROL_ONLY = {"Sm_core"}
IMMUNE_PANELS = [p for p in PANELS if p not in CONTROL_ONLY]
ALL_GENES = sorted({g for v in PANELS.values() for g in v})


def build_map():
    done = {}
    if os.path.exists(MAPF):
        m = pd.read_csv(MAPF, sep="\t", header=None, names=["sym", "ensg"], dtype=str)
        m["ensg"] = m["ensg"].fillna("")
        for _, r in m.iterrows():
            if r["ensg"] and r["sym"] not in done:
                done[r["sym"]] = r["ensg"]
    todo = [g for g in ALL_GENES if not done.get(g)]
    if todo:
        op = urllib.request.build_opener(urllib.request.ProxyHandler({}))

        def q(g):
            u = ("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"
                 f"{urllib.parse.quote(g)}?content-type=application/json")
            for _ in range(4):
                try:
                    with op.open(u, timeout=30) as r:
                        v = pd.read_json(r).to_dict("records")
                    c = sorted({str(x.get("id", "")) for x in v
                                if str(x.get("id", "")).startswith("ENSG")})
                    return g, ",".join(c)
                except Exception:
                    time.sleep(2)
            return g, ""

        with ThreadPoolExecutor(6) as ex:
            for g, e in ex.map(q, todo):
                done[g] = e
                with open(MAPF, "a") as fh:
                    fh.write(f"{g}\t{e}\n")
    miss = [g for g in ALL_GENES if not done.get(g)]
    if miss:
        print(f"  [warn] 未解析: {miss}")
    return done


def scan_expr(ensg_set):
    """一次扫描 740MB：目标基因行 + 每样本全局表达均值。"""
    rows, t0 = {}, time.time()
    with gzip.open(EXPR, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = np.array(header[1:])
        n = len(cols)
        s_sum = np.zeros(n, dtype=np.float64)
        n_gene = 0
        for line in f:
            gid, rest = line.split("\t", 1)
            gid = gid.split(".")[0]
            try:
                vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            except ValueError:
                continue
            if len(vals) != n:
                continue
            n_gene += 1
            s_sum += vals
            if gid in ensg_set:
                rows[gid] = vals
    print(f"  [scan] {n_gene} genes x {n} samples, {time.time()-t0:.0f}s; "
          f"命中 {len(rows)}/{len(ensg_set)}", flush=True)
    stat = pd.DataFrame({"TotalExpr": s_sum / max(n_gene, 1)}, index=cols)
    return pd.DataFrame(rows, index=cols), stat


def rank_partial(x, y, covs):
    """秩偏相关：rank -> 对协变量 OLS -> 残差 Pearson。k=0 时 = Spearman。"""
    n = len(x)
    xr = st.rankdata(x)
    yr = st.rankdata(y)
    k = covs.shape[1] if covs is not None and covs.size else 0
    if k == 0:
        r = float(np.corrcoef(xr, yr)[0, 1])
    else:
        X = np.column_stack([np.ones(n), st.rankdata(covs, axis=0)])
        rx = xr - X @ np.linalg.lstsq(X, xr, rcond=None)[0]
        ry = yr - X @ np.linalg.lstsq(X, yr, rcond=None)[0]
        r = float(np.corrcoef(rx, ry)[0, 1])
    df = n - 2 - k
    if df <= 0 or not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan, df
    t = r * np.sqrt(df / max(1 - r * r, 1e-12))
    return r, float(2 * st.t.sf(abs(t), df)), df


def fisher_ci(r, n):
    if not np.isfinite(r) or abs(r) >= 1 or n < 5:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    gmap = build_map()
    ensg = {x for g in ALL_GENES for x in gmap.get(g, "").split(",") if x}
    E, stat = scan_expr(ensg)

    # 样本 key：TCGA-xx-xxxx-01 -> TCGA-xx-xxxx
    # 同一病人可能有多个 aliquot（-01 原发 / -06 转移 / -11 正常），
    # 必须先按「原发优先」去重，否则 join 会膨胀（8928 > 8260）。
    def to_patient(df):
        full = pd.Index(df.index)
        key = pd.Index([s[:12] for s in full])
        pri = pd.Series([0 if str(s).endswith("-01") else 1 for s in full],
                        index=range(len(full)))
        keep = (pd.DataFrame({"k": key, "p": pri.values})
                .sort_values("p").groupby("k").head(1).index)
        out = df.iloc[keep].copy()
        out.index = key[keep]
        return out

    E = to_patient(E)
    stat = to_patient(stat)

    bd = pd.read_csv(BURDEN, sep="\t")
    bd["patient"] = bd["patient"].astype(str).str[:12]
    bd = bd.drop_duplicates("patient")
    bd["cohort"] = bd["projectID"].str.replace("TCGA-", "", regex=False)
    print(f"  burden: {len(bd)} 样本 / {bd.cohort.nunique()} 癌种")

    # 面板分数
    scores = {}
    for pname, genes in PANELS.items():
        es = [x for g in genes for x in gmap.get(g, "").split(",") if x]
        es = [x for x in es if x in E.columns]
        scores[pname] = E[es].mean(axis=1)
        print(f"  {pname:<18} {len(es)} 基因")
    S = pd.DataFrame(scores)
    S["TotalExpr"] = stat["TotalExpr"].reindex(S.index)
    S = S[~S.index.duplicated()]
    S = S.join(bd.set_index("patient")[["RI_burden", "cohort"]], how="inner")
    S = S.dropna(subset=["RI_burden"])
    print(f"  合并后样本: {len(S)}  癌种: {S.cohort.nunique()}")

    # 综合免疫浸润分数 = 免疫面板 z-score 的均值（癌种内 z）
    Z = S.groupby("cohort")[IMMUNE_PANELS].transform(
        lambda v: (v - v.mean()) / v.std(ddof=0))
    S["Immune_aggregate"] = Z.mean(axis=1)
    PANELS_OUT = IMMUNE_PANELS + ["Immune_aggregate", "Sm_core"]

    rows = []
    samp = []
    for coh, g in S.groupby("cohort"):
        if len(g) < MIN_N:
            print(f"  [排除] {coh} n={len(g)}")
            continue
        for pn in PANELS_OUT:
            y = g["RI_burden"].values
            x = g[pn].values
            r0, p0, _ = rank_partial(x, y, None)
            r1, p1, df = rank_partial(x, y, g[["TotalExpr"]].values)
            lo0, hi0 = fisher_ci(r0, len(g))
            lo1, hi1 = fisher_ci(r1, len(g) - 1)
            rows.append(dict(cohort=coh, panel=pn, n=len(g),
                             rho_unadj=r0, ci_lo=lo0, ci_hi=hi0, p_unadj=p0,
                             rho_adj=r1, ci_lo_adj=lo1, ci_hi_adj=hi1, p_adj=p1))
        gg = g.assign(**{f"z_{p}": Z.loc[g.index, p] for p in IMMUNE_PANELS})
        samp.append(gg[["cohort", "RI_burden", "TotalExpr", "Immune_aggregate"]
                       + IMMUNE_PANELS + ["Sm_core"]])

    R = pd.DataFrame(rows)
    # 面板内跨癌种 BH
    R["q_unadj"] = np.nan
    R["q_adj"] = np.nan
    for pn, g in R.groupby("panel"):
        R.loc[g.index, "q_unadj"] = bh(g["p_unadj"].values)
        R.loc[g.index, "q_adj"] = bh(g["p_adj"].values)
    R.to_csv(os.path.join(OUT, "pancan_tme_map.tsv"), sep="\t", index=False)
    pd.concat(samp).to_csv(os.path.join(OUT, "pancan_tme_samples.tsv"),
                           sep="\t", index=False)
    print("\n-> pancan_tme_map.tsv / pancan_tme_samples.tsv")

    print("\n=== 综合免疫浸润 vs IR burden（校正后 rho_adj，按升序）===")
    t = R[R.panel == "Immune_aggregate"].sort_values("rho_adj")
    print(t[["cohort", "n", "rho_unadj", "rho_adj", "q_adj"]].round(3).to_string(index=False))
    print(f"\n负相关癌种数: {(t.rho_adj < 0).sum()}/{len(t)}; "
          f"q<0.05: {(t.q_adj < 0.05).sum()}/{len(t)}")
    print("\n=== Sm_core 对照（校正后）===")
    s = R[R.panel == "Sm_core"].sort_values("rho_adj")
    print(s[["cohort", "n", "rho_unadj", "rho_adj", "q_adj"]].round(3).to_string(index=False))
    print("\n=== 各面板跨癌种中位数 ===")
    print(R.groupby("panel")[["rho_unadj", "rho_adj"]].median().round(3).to_string())


if __name__ == "__main__":
    main()
