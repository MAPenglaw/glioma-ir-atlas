#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""45 — TME 耦合的混杂排查：偏相关 + 总表达量校正。

为什么必须有这一步
------------------
1. 旧口径（OncoSplicing, n=151）里 GBM 的「全基因平均表达 MeanExprAll」与
   RI burden 相关 rho=0.352 (p=9e-6)，是当时最强的混杂；但现行 Kahles 口径
   （GBM 150 / LGG 498）**从未检验过这个协变量**。若不校正就宣称
   「LGG 是免疫冷表型」，审稿人一句「这只是总表达量的影子」就能推翻。
2. 44 号发现的 LGG 11/16 显著负相关，必须证明它不是
   (a) 总表达量 / 检出深度的技术混杂，
   (b) 正常脑组织（神经元）污染的组成效应，
   (c) IDH 分层的替身。
3. 顺带对主线结论做压力测试：剪接体 Sm core 与 burden 的关联，
   在校正总表达量后是否仍然站得住。

方法
----
- 协变量（癌种内 z-score / 秩正态）：
    TotalExpr      = 样本全部基因 log2(TPM+0.001) 的均值（技术混杂主代理）
    TotalExpr_det  = 检出基因比例（log2TPM > 0 的基因占比）
    Neuronal       = 正常脑/神经元污染
    Sm_core / Proliferation / Tumor_content = 44 号的既有对照
    IDH_status     = IDHmut / IDHwt（分型替身，来自 33 号 cBioPortal）
- 偏相关 = 对 x、y 各自做「秩正态(inverse-normal) → 对协变量 OLS → 取残差」，
  再对残差做 Pearson；df = n - 2 - k。这是 Spearman 型偏相关的标准实现。
- 模型梯度：M0(未校正) → M1(+TotalExpr) → M2(+Neuronal) → M3(+剪接体/增殖)
  → M4(+纯度) → M5(+IDH)。看信号衰减轨迹，而不是只看一个数字。

输出
----
  output_pancan/tme_partial_coupling.tsv   癌种 × 面板 × 模型 的偏 rho/CI/p
  output_pancan/tme_confounders.tsv        样本级协变量（供绘图）
  output_pancan/tme_partial_samples.tsv    样本级 面板分数 + burden（供绘图）
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
MOL = os.path.join(OUT, "mol_burden_plot.tsv")

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
    "Neuronal_normal_brain": ["RBFOX3", "SYT1", "SNAP25", "GRIN1", "GABRA1",
                              "MAP2", "NEFL", "NEFM", "STMN2", "SYP"],
    "TotalExpr": [],      # 由全基因扫描得到，不是基因面板
    "TotalExpr_det": [],
}
GROUPS = {**TME_GROUPS, **CONTROL_GROUPS}
GENE_GROUPS = {k: v for k, v in GROUPS.items() if v}
ALL_GENES = sorted({g for v in GENE_GROUPS.values() for g in v})

# 面板覆盖基因太少时给出提示（<3 个基因的面板不进主结论）
MIN_GENES = 3


def build_map():
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
    with ThreadPoolExecutor(12) as ex:
        for i, (g, v) in enumerate(ex.map(q, todo), 1):
            done[g] = v
            if v:
                with open(MAPF, "a") as fh:
                    fh.write(f"{g}\t{v}\n")
            if i % 20 == 0:
                print(f"  [map] {i}/{len(todo)}", flush=True)
    return done


def scan_expr(ensg_set, sample_frac=0.05, max_sample_genes=6000):
    """一次扫描 740MB：同时拿目标基因行 + 每样本全局表达统计。

    TotalExpr      = 全部基因 log2(TPM+0.001) 的均值（技术混杂主代理）
    TotalExpr_det  = 检出基因比例（log2TPM > 0 的基因占比）
    TotalExpr_hi/lo= **按基因表达水平分层**后的样本均值。
        用途：判别「burden ~ TotalExpr」到底是技术还是生物学。
        若相关主要由低表达基因层驱动 → 检出/深度膨胀（技术）；
        若两层贡献相当 → 整体表达偏移（更可能是组成/生物学）。
    """
    rows, t0 = {}, time.time()
    rng = np.random.default_rng(20260928)
    with gzip.open(EXPR, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = np.array(header[1:])
        n = len(cols)
        s_sum = np.zeros(n, dtype=np.float64)
        s_det = np.zeros(n, dtype=np.float64)
        n_gene = 0
        pool = []          # 抽样的基因行（用于事后分层）
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
            s_det += (vals > 0)
            if gid in ensg_set:
                rows[gid] = vals
            if len(pool) < max_sample_genes and rng.random() < sample_frac:
                pool.append(vals.copy())
    print(f"  [scan] {n_gene} genes x {n} samples, {time.time()-t0:.0f}s; "
          f"target rows {len(rows)}/{len(ensg_set)}; pooled {len(pool)}",
          flush=True)
    stat = pd.DataFrame({"TotalExpr": s_sum / max(n_gene, 1),
                         "TotalExpr_det": s_det / max(n_gene, 1)},
                        index=cols)
    if pool:
        Pm = np.vstack(pool)                 # genes x samples
        gmean = Pm.mean(axis=1)
        cut = np.median(gmean)
        hi = gmean > cut
        stat["TotalExpr_hi"] = Pm[hi].mean(axis=0)
        stat["TotalExpr_lo"] = Pm[~hi].mean(axis=0)
        print(f"  [scan] 分层: hi n={int(hi.sum())} (mean log2TPM "
              f"{gmean[hi].mean():.2f}), lo n={int((~hi).sum())} "
              f"({gmean[~hi].mean():.2f})", flush=True)
    return pd.DataFrame(rows, index=cols), stat


def rankit(x):
    """秩变换。对秩做 Pearson == Spearman，故 M0 与 44 号口径完全一致。"""
    return st.rankdata(x).astype(float)


def invnorm(x):
    r = st.rankdata(x)
    return st.norm.ppf(r / (len(r) + 1.0))


def partial_corr(x, y, covs):
    """对协变量 OLS → 残差 Pearson。x/y/covs 需已按所选尺度变换。"""
    n = len(x)
    k = covs.shape[1] if covs is not None and covs.size else 0
    if k == 0:
        r = float(np.corrcoef(x, y)[0, 1])
    else:
        X = np.column_stack([np.ones(n), covs])
        bx = np.linalg.lstsq(X, x, rcond=None)[0]
        by = np.linalg.lstsq(X, y, rcond=None)[0]
        rx = x - X @ bx
        ry = y - X @ by
        r = float(np.corrcoef(rx, ry)[0, 1])
    df = n - 2 - k
    if df <= 0 or not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan, df
    t = r * np.sqrt(df / max(1 - r * r, 1e-12))
    p = 2 * st.t.sf(abs(t), df)
    return r, float(p), df


def fisher_ci(r, n):
    if n < 10 or not np.isfinite(r) or abs(r) >= 0.999:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(max(n - 3 - 1, 1))
    return (np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se))


def mde(n, k=0, alpha=0.05, power=0.8):
    za, zb = st.norm.ppf(1 - alpha / 2), st.norm.ppf(power)
    return float(np.tanh((za + zb) / np.sqrt(max(n - 3 - k, 1))))


# 注意：TotalExpr 与 TotalExpr_det 的 Spearman 相关高达 0.91–0.92，
# 二者**不可同时入模**（共线会让偏相关系数失控地跳变）。故改为互为替代。
MODELS = [
    ("M0_raw", []),
    ("M1_TotalExpr", ["TotalExpr"]),
    ("M2_plus_Neuronal", ["TotalExpr", "Neuronal_normal_brain"]),
    ("M3_plus_Spliceosome", ["TotalExpr", "Neuronal_normal_brain",
                             "Sm_core", "Proliferation"]),
    ("M4_plus_Purity", ["TotalExpr", "Neuronal_normal_brain", "Sm_core",
                        "Proliferation", "Tumor_content"]),
    ("M5_plus_IDH", ["TotalExpr", "Neuronal_normal_brain", "Sm_core",
                     "Proliferation", "Tumor_content", "IDHmut"]),
    # 敏感性 A：换用「检出广度」替代「平均表达」作技术代理
    ("M6_alt_Detection", ["TotalExpr_det", "Neuronal_normal_brain", "Sm_core",
                          "Proliferation", "Tumor_content", "IDHmut"]),
    # 敏感性 B：把技术代理拆成高/低表达两层同时入模（最保守的技术校正）
    ("M7_hi_lo_layers", ["TotalExpr_hi", "TotalExpr_lo",
                         "Neuronal_normal_brain", "Sm_core", "Proliferation",
                         "Tumor_content", "IDHmut"]),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    gmap = build_map()
    cands = {g: [x for x in gmap.get(g, "").split(",") if x] for g in ALL_GENES}
    ensg_set = {x for v in cands.values() for x in v}
    expr, glob = scan_expr(ensg_set)

    chosen, used = {}, {}
    for g, cs in cands.items():
        ok = [c for c in cs if c in expr.columns]
        if not ok:
            chosen[g] = None
            continue
        chosen[g] = max(ok, key=lambda c: float(expr[c].mean()))
    keep_cols = {v for v in chosen.values() if v}
    expr = expr[[c for c in expr.columns if c in keep_cols]]
    inv = {}
    for g, v in chosen.items():
        if v and v not in inv:
            inv[v] = g
    expr.columns = [inv[c] for c in expr.columns]
    for gname, genes in GENE_GROUPS.items():
        used[gname] = [g for g in genes if g in expr.columns]
    thin = {k: v for k, v in used.items() if len(v) < MIN_GENES}
    if thin:
        print(f"  [warn] 基因数 <{MIN_GENES} 的面板: {thin}")

    # 只保留原发肿瘤(-01)，折叠到 patient key
    keep = [i for i, s in enumerate(expr.index) if s.endswith("-01")]
    expr = expr.iloc[keep]
    expr.index = [s[:12] for s in expr.index]
    expr = expr[~expr.index.duplicated(keep="first")]
    glob = glob.iloc[keep]
    glob.index = [s[:12] for s in glob.index]
    glob = glob[~glob.index.duplicated(keep="first")]

    mol = pd.read_csv(MOL, sep="\t")
    mol["pt"] = mol["sample"].astype(str).str[:12]
    mol = mol.drop_duplicates("pt").set_index("pt")

    rows, samp_rows = [], []
    for cohort, bpath in BURDEN.items():
        bd = pd.read_csv(bpath, sep="\t", index_col=0)
        bd.index = [s[:12] for s in bd.index]
        bd = bd[~bd.index.duplicated(keep="first")]
        common = bd.index.intersection(expr.index).intersection(glob.index)
        print(f"\n=== {cohort}: burden n={len(bd)}, matched n={len(common)} ===")
        e = expr.loc[common]
        g = glob.loc[common]
        y_raw = bd.loc[common, "RI_burden"].values.astype(float)

        # ---- 样本级分数矩阵 ----
        S = pd.DataFrame(index=common)
        for gname, genes in GENE_GROUPS.items():
            es = used[gname]
            if len(es) < 2:
                continue
            sub = e[es].astype(float)
            z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
            S[gname] = z.mean(axis=1).values
        # 存原始值；尺度变换（秩 / 秩正态）留到分析时统一做
        S["TotalExpr"] = g["TotalExpr"].values
        S["TotalExpr_det"] = g["TotalExpr_det"].values
        for c in ("TotalExpr_hi", "TotalExpr_lo"):
            if c in g.columns:
                S[c] = g[c].values
        S["RI_burden"] = y_raw
        ok = S.notna().all(axis=1)
        S = S[ok]
        print(f"  完整样本 n={len(S)}")

        # IDH（缺失填 0，并在结果里记下覆盖率）
        idh = mol["IDH_status"].reindex(S.index)
        idhmut = (idh == "IDHmut").astype(float)
        cov_idh = idh.notna().astype(float).mean()
        S["IDHmut"] = idhmut.fillna(0.0).values
        print(f"  IDH 覆盖: {cov_idh:.1%}  ({int((idh=='IDHmut').sum())} mut / "
              f"{int((idh=='IDHwt').sum())} wt)")

        yr = rankit(S["RI_burden"].values)      # 秩尺度（主口径）
        yi = invnorm(S["RI_burden"].values)     # 秩正态尺度（敏感性）
        # ---- 未校正的 burden ~ 总表达量（关键混杂量级）----
        r_te, p_te = st.spearmanr(S["TotalExpr"], S["RI_burden"])
        print(f"  [关键] RI_burden ~ TotalExpr: rho={r_te:+.3f} p={p_te:.2e}")
        for c in ("TotalExpr_hi", "TotalExpr_lo"):
            if c in S.columns:
                rh, ph = st.spearmanr(S[c], S["RI_burden"])
                print(f"         RI_burden ~ {c}: rho={rh:+.3f} p={ph:.2e}")

        for gname in [c for c in S.columns
                      if c not in ("RI_burden", "IDHmut")]:
            rs, ps = st.spearmanr(S[gname].values, S["RI_burden"].values)
            xr = rankit(S[gname].values.astype(float))
            xi = invnorm(S[gname].values.astype(float))
            for mname, covs in MODELS:
                cs = [c for c in covs if c != gname and c in S.columns]
                if cs:
                    Craw = S[cs].values.astype(float)
                    Cr = np.column_stack([rankit(Craw[:, j])
                                          for j in range(Craw.shape[1])])
                    Ci = np.column_stack([invnorm(Craw[:, j])
                                          for j in range(Craw.shape[1])])
                else:
                    Cr = Ci = np.empty((len(S), 0))
                r, p, df = partial_corr(xr, yr, Cr)
                ri, pi, _ = partial_corr(xi, yi, Ci)
                if mname == "M0_raw":       # M0 直接用 Spearman 的精确 p
                    r, p = float(rs), float(ps)
                lo, hi = fisher_ci(r, len(S) - len(cs))
                rows.append(dict(cohort=cohort, panel=gname, model=mname,
                                 rho=r, ci_lo=lo, ci_hi=hi, p=p,
                                 rho_inv=ri, p_inv=pi,
                                 rho_raw_spearman=float(rs),
                                 p_raw_spearman=float(ps),
                                 n=len(S), k_cov=len(cs),
                                 mde80=mde(len(S), len(cs))))
            samp_rows.append(S.assign(cohort=cohort, panel=gname))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "tme_partial_coupling.tsv"),
               sep="\t", index=False)
    pd.concat([s.reset_index().rename(columns={"index": "sample"})
               for s in samp_rows]) \
        .to_csv(os.path.join(OUT, "tme_partial_samples.tsv"),
                sep="\t", index=False)
    S_all = pd.concat([s.reset_index().rename(columns={"index": "sample"})
                       for s in samp_rows])
    ccols = ["sample", "cohort", "RI_burden", "TotalExpr", "TotalExpr_det",
             "TotalExpr_hi", "TotalExpr_lo", "Neuronal_normal_brain",
             "Sm_core", "Proliferation", "Tumor_content", "IDHmut"]
    conf = S_all.drop_duplicates("sample")
    conf = conf[[c for c in ccols if c in conf.columns]]
    conf.to_csv(os.path.join(OUT, "tme_confounders.tsv"), sep="\t", index=False)

    print("\n=== 关键面板：模型梯度下的偏 rho（* p<0.05, ** p<0.001）===")
    key = ["TotalExpr", "TotalExpr_hi", "TotalExpr_lo", "Sm_core", "CD8_T",
           "IFN_gamma_Ayers", "Endothelial", "Cytotoxic", "Proliferation"]
    for coh in ["GBM", "LGG"]:
        sub = res[(res.cohort == coh) & (res.panel.isin(key))].copy()
        if sub.empty:
            continue
        piv = sub.pivot(index="panel", columns="model", values="rho")
        ppv = sub.pivot(index="panel", columns="model", values="p")
        piv = piv[[m for m, _ in MODELS if m in piv.columns]]
        ppv = ppv[[m for m, _ in MODELS if m in ppv.columns]]
        lab = piv.copy().astype(object)
        for c in piv.columns:
            lab[c] = [f"{v:+.3f}{'**' if p < 0.001 else ('*' if p < 0.05 else '')}"
                      for v, p in zip(piv[c], ppv[c])]
        print(f"\n--- {coh} ---")
        print(lab.to_string())


if __name__ == "__main__":
    main()
