#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""20 — 把 G3（IR 事件基因组序列特征）迁移到 Kahles 2018 事件集。

与旧版 11 号脚本的差异：
  - 数据源：KM_IR_GBM.psi.gz（Kahles 2018，5 列 meta + 样本 PSI），事件全集 123,390
  - 坐标格式：event_coordinates = "e1s:e1e:e2s:e2e"，alt_region = "is:ie"（1-based inclusive）
  - 无 strand 字段 → 用内含子两端序列自证（GT..AG = + 链；CT..AC = - 链）
  - 无预计算 MedianPSI/stdPSI → 从 GBM 样本 PSI 自算
  - 无 isoform 注释 → 去掉 n_iso 特征

三类（Braunschweig 2014）：Variable = stdPSI≥0.10；其余按 MedianPSI 三分位取 Spliced/Retained。
"""
import os, json, time, urllib.request
import numpy as np
import pandas as pd

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
DATA = os.path.join(ROOT, "ri_atlas", "data")
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "g3_km_intron_seqs.json")
API = "https://grch37.rest.ensembl.org/sequence/region/human?content-type=application/json"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s):
    return s.translate(COMP)[::-1]

np.random.seed(20260927)

# ---------------------------------------------------------------- 1. 读 PSI + 坐标
print("[1] 读 Kahles GBM IR PSI，解析坐标 ...", flush=True)
df = pd.read_csv(os.path.join(DATA, "KM_IR_GBM.psi.gz"), sep="\t",
                 compression="gzip", low_memory=False)
meta = df[["event_id", "event_type", "event_chr", "event_coordinates", "alt_region_coordinates"]]
vals = df.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
print(f"    事件 {len(df):,} × 样本 {vals.shape[1]}")

# 解析坐标
ec = meta["event_coordinates"].str.split(":", expand=True).astype(float)
alt = meta["alt_region_coordinates"].str.split(":", expand=True).astype(float)
meta["e1s"], meta["e1e"] = ec[0], ec[1]
meta["e2s"], meta["e2e"] = ec[2], ec[3]
meta["i_s"], meta["i_e"] = alt[0], alt[1]
meta["chrom"] = meta["event_chr"].map(str).str.replace("chr", "", regex=False)
meta["ilen"] = meta["i_e"] - meta["i_s"] + 1
meta["e1len"] = meta["e1e"] - meta["e1s"] + 1
meta["e2len"] = meta["e2e"] - meta["e2s"] + 1

# 自算 MedianPSI / stdPSI / 检出率
meta["MedianPSI"] = vals.median(axis=1)
meta["stdPSI"] = vals.std(axis=1)
meta["pctSamples"] = vals.notna().mean(axis=1)

# 过滤：检出率足够（能可靠算 median/std）
MIN_PCT = 0.5
use = meta[meta["pctSamples"] >= MIN_PCT].copy()
print(f"    检出率≥{MIN_PCT} 的事件：{len(use):,}")

# ---------------------------------------------------------------- 2. 三组分类
VAR_CUT = 0.10
grp = pd.Series("Intermediate", index=use.index)
grp[use["stdPSI"] >= VAR_CUT] = "Variable"
rest = use[use["stdPSI"] < VAR_CUT]
p = rest["MedianPSI"].quantile([1 / 3, 2 / 3]).values
grp[(use["stdPSI"] < VAR_CUT) & (use["MedianPSI"] <= p[0])] = "Spliced"
grp[(use["stdPSI"] < VAR_CUT) & (use["MedianPSI"] >= p[1])] = "Retained"
use["grp"] = grp
print("    三组构成：\n" + use["grp"].value_counts().to_string())
print(f"    MedianPSI 三分位切点：{p[0]:.3f} / {p[1]:.3f}")

# 只对三组（非 Intermediate）取序列
target = use[use["grp"].isin(["Spliced", "Retained", "Variable"])].copy()
print(f"    需取序列的三组事件：{len(target):,}")

# ---------------------------------------------------------------- 3. 取序列
# 禁用代理直连 Ensembl（本机代理对 grch37.rest.ensembl.org 间歇性 404/502）
_NO_PROXY_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def fetch_batch(regions, tries=4):
    req = urllib.request.Request(API, data=json.dumps({"regions": regions}).encode(),
                                 headers={"Content-Type": "application/json"})
    for t in range(tries):
        try:
            with _NO_PROXY_OPENER.open(req, timeout=180) as r:
                return json.loads(r.read())
        except Exception as e:
            if t == tries - 1:
                body = ""
                try:
                    body = e.read().decode()[:200] if hasattr(e, "read") else str(e)
                except Exception:
                    body = str(e)
                print(f"    ! 批次失败 {e} | body: {body}", flush=True)
                return None
            time.sleep(2 + 3 * t)
    return None


def run_fetch(target):
    seqs = {}
    if os.path.exists(CACHE):
        seqs = json.load(open(CACHE))
        print(f"    [缓存] 已有 {len(seqs):,} 条，续传", flush=True)
    seqs, keys, regs = seqs, [], []
    for k, r in target.iterrows():
        if k in seqs:
            continue
        s0 = max(0, int(r["i_s"]) - 3)
        e0 = int(r["i_e"]) + 3
        keys.append(k)
        regs.append(f"{r['chrom']}:{s0}-{e0}")   # Ensembl 1-based inclusive
    print(f"    待取 {len(regs):,} 个（已缓存 {len(seqs):,}）", flush=True)
    BATCH = 45
    batches = [list(range(i, min(i + BATCH, len(regs)))) for i in range(0, len(regs), BATCH)]
    for bi, bidx in enumerate(batches):
        res = fetch_batch([regs[j] for j in bidx])
        if res is None:
            continue
        for j, item in zip(bidx, res):
            seqs[keys[j]] = item.get("seq", "")
        if (bi + 1) % 100 == 0 or bi == len(batches) - 1:
            json.dump(seqs, open(CACHE, "w"))
            print(f"      批 {bi+1}/{len(batches)} 已获取 {len(seqs):,}（已存缓存）", flush=True)
        time.sleep(0.05)
    json.dump(seqs, open(CACHE, "w"))
    return seqs

print("[2] Ensembl GRCh37 取序列（内含子 ±3bp）...", flush=True)
seqs = run_fetch(target)

# ---------------------------------------------------------------- 4. 特征（含 strand 自证）
print("[3] 计算特征 + strand 自证 ...", flush=True)
rows, drop = [], 0
for k, r in target.iterrows():
    raw = seqs.get(k, "")
    if not raw:
        drop += 1
        continue
    ilen = int(r["ilen"])
    raw = raw.upper()
    if len(raw) != ilen + 6:
        drop += 1
        continue
    # strand 自证：+ 链内含子两端 GT..AG；- 链为 CT..AC（互补）
    ss5_raw, ss3_raw = raw[3:5], raw[3 + ilen - 2:3 + ilen]
    if ss5_raw == "CT" and ss3_raw == "AC":
        raw = revcomp(raw)
    intron = raw[3:3 + ilen]
    if len(intron) < 20:
        drop += 1
        continue
    gc = (intron.count("G") + intron.count("C")) / len(intron)
    ss5, ss3 = intron[:2], intron[-2:]
    pt20 = intron[-20:]
    y20 = (pt20.count("C") + pt20.count("T")) / 20
    tail = intron[-25:]
    maxy, cur = 0, 0
    for ch in tail:
        cur = cur + 1 if ch in "CT" else 0
        maxy = max(maxy, cur)
    rows.append({
        "event": k, "chrom": r["chrom"], "ilen": ilen,
        "e1len": int(r["e1len"]), "e2len": int(r["e2len"]),
        "gc": gc, "ss5": ss5, "ss3": ss3,
        "ss5_class": "GT" if ss5 == "GT" else ("GC" if ss5 == "GC" else ("AT" if ss5 == "AT" else "other")),
        "y20": y20, "maxY": maxy,
        "MedianPSI": r["MedianPSI"], "stdPSI": r["stdPSI"], "pctSamples": r["pctSamples"],
        "grp": r["grp"],
    })
F = pd.DataFrame(rows)
print(f"    可用于分析 {len(F):,}｜丢弃 {drop:,}")
F.to_csv(os.path.join(OUT, "G3_km_intron_features.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- 5. 统计
from scipy import stats
NUM = ["ilen", "gc", "e1len", "e2len", "y20", "maxY"]
print("\n[5] 三类特征对比（中位数，KW p，S vs R p）")
print(f"{'feature':10s} {'Spliced':>12s} {'Retained':>12s} {'Variable':>12s} {'KW p':>12s} {'SvsR p':>12s}")
for c in NUM:
    g = [F.loc[F["grp"] == k, c].dropna() for k in ["Spliced", "Retained", "Variable"]]
    if min(len(x) for x in g) < 10:
        continue
    kw = stats.kruskal(*g).pvalue
    sr = stats.mannwhitneyu(g[0], g[1], alternative="two-sided").pvalue
    print(f"{c:10s} {g[0].median():12.1f} {g[1].median():12.1f} {g[2].median():12.1f} {kw:12.2e} {sr:12.2e}")

print("\n[6] 5' 剪接位点类别构成（%）")
ct = pd.crosstab(F["grp"], F["ss5_class"], normalize="index").mul(100).round(2)
print(ct.reindex(["Spliced", "Retained", "Variable"]).to_string())
sub = F[F["grp"].isin(["Spliced", "Retained"])]
tab = pd.crosstab(sub["grp"], sub["ss5_class"])[["GT", "GC", "AT", "other"]]
chi2 = stats.chi2_contingency(tab.values)
print(f"  非经典 5'ss Spliced vs Retained：卡方 p = {chi2[1]:.3e}")
print(f"  GC-AG 占比：Spliced {ct.loc['Spliced','GC']:.2f}% vs Retained {ct.loc['Retained','GC']:.2f}%")
gc_fold = ct.loc["Retained", "GC"] / max(ct.loc["Spliced", "GC"], 1e-9)
print(f"  GC-AG 富集倍数：{gc_fold:.1f}×")

# 长度匹配条件分析
print("\n[7] 长度匹配（十分位分层）后的条件检验")
F["dlen"] = pd.qcut(F["ilen"], 10, labels=False, duplicates="drop")
for c in ["gc", "y20", "maxY", "e1len", "e2len"]:
    S, R = F[F["grp"] == "Spliced"], F[F["grp"] == "Retained"]
    ds, ps = [], []
    for d in sorted(set(F["dlen"].dropna())):
        a = S.loc[S["dlen"] == d, c].dropna(); b = R.loc[R["dlen"] == d, c].dropna()
        if len(a) < 30 or len(b) < 30:
            continue
        ps.append(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    if ps:
        # 加权 Z（Stouffer）近似
        import scipy.stats as st
        zs = st.norm.ppf(1 - np.array(ps))
        z_comb = zs.sum() / np.sqrt(len(zs))
        p_comb = 2 * (1 - st.norm.cdf(abs(z_comb)))
        print(f"  {c:6s}: 各十分位 p 中位 {np.median(ps):.2e}，合并 p = {p_comb:.2e}")
print("\nDONE")
