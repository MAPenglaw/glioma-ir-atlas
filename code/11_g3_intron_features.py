#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G3：IR 事件的基因组序列特征（"atlas" 标配）

设计（依据 Braunschweig 2014 GR 对 IR 的分类框架）：
  按「滞留倾向 × 变异度」把 IR 事件分为三类：
    - Retained (构成性滞留)：MedianPSI 高、stdPSI 低
    - Spliced  (高效剪除)  ：MedianPSI 低、stdPSI 低
    - Variable (调控性滞留)：stdPSI 高
  比较三类的内含子长度 / GC / 剪接位点类别 / 多嘧啶区 / 侧翼外显子长度 / isoform 复杂度

数据源：
  - SplAdder_TCGA_GBM_PSI.txt.gz  → EventRegion 坐标（0-based, half-open）
  - Ensembl REST GRCh37 在线取序列（本地无基因组 fasta）

版本自证：+ 链内含子两端必须是 GT..AG —— 已验证 GRCh37 命中（GRCh38 全错）。
"""
import os, json, sys, time, urllib.request, urllib.error
import numpy as np
import pandas as pd

np.random.seed(20260921)

SRC = ("/Users/peng/课题/SNRPC-UQCRC1_2027申报储备_v1/"
       "07_原始项目数据_2026_splicing_analysis/2026 splicing analysis/splicing_analysis")
OUT = os.path.expanduser("/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output")
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "g3_intron_seqs.json")
API = "https://grch37.rest.ensembl.org/sequence/region/human?content-type=application/json"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- 1. 坐标
print("[1] 解析 SplAdder IR 事件坐标 ...", flush=True)
meta = pd.read_csv(f"{SRC}/SplAdder_TCGA_GBM_PSI.txt.gz", sep="\t", index_col=0,
                   compression="gzip", low_memory=False, nrows=300000)
ir = meta[meta["SpliceType"] == "IR"].copy()

parts = ir["EventRegion"].str.split("-", expand=True)
def_seg = lambda col, i: col.str.split(":", expand=True)[i].astype(float)
ir["e1s"], ir["e1e"] = def_seg(parts[0], 0), def_seg(parts[0], 1)
ir["i_s"], ir["i_e"] = def_seg(parts[1], 0), def_seg(parts[1], 1)
ir["e2s"], ir["e2e"] = def_seg(parts[2], 0), def_seg(parts[2], 1)
ir["chrom"] = ir["ChrStrand"].str.partition("|")[0].str.strip()
ir["strand"] = ir["ChrStrand"].str.partition("|")[2].str.strip()
ir["ilen"] = ir["i_e"] - ir["i_s"]
ir["e1len"] = ir["e1e"] - ir["e1s"]
ir["e2len"] = ir["e2e"] - ir["e2s"]
ir = ir.dropna(subset=["i_s", "i_e", "ilen"])
print(f"    IR 事件 {len(ir):,}｜长度中位数 {ir['ilen'].median():.0f} bp")


# ---------------------------------------------------------------- 2. 取序列
def fetch_batch(regions, tries=4):
    req = urllib.request.Request(API, data=json.dumps({"regions": regions}).encode(),
                                 headers={"Content-Type": "application/json"})
    for t in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read())
        except Exception as e:
            if t == tries - 1:
                print(f"    ! 批次失败 {e}", flush=True)
                return None
            time.sleep(2 + 3 * t)
    return None


def run_fetch(ir):
    """每次请求 45 个 region；返回 {event_id: seq(+内含子,转录本方向)}"""
    if os.path.exists(CACHE):
        seqs = json.load(open(CACHE))
        print(f"    [缓存] 已存在 {len(seqs):,} 条序列，跳过下载")
        return seqs
    seqs, keys, regs = {}, [], []
    for k, r in ir.iterrows():
        # 请求区间 = 上游外显子末3bp + 内含子 + 下游外显子首3bp（0-based）
        s0 = max(0, int(r["i_s"]) - 3)
        e0 = int(r["i_e"]) + 3
        keys.append(k)
        regs.append(f"{r['chrom'].replace('chr','')}:{s0+1}-{e0}")   # Ensembl 1-based inclusive
    batches = [list(range(i, min(i + 45, len(regs)))) for i in range(0, len(regs), 45)]
    print(f"    {len(regs):,} 个区域 → {len(batches)} 批", flush=True)
    for bi, bidx in enumerate(batches):
        res = fetch_batch([regs[j] for j in bidx])
        if res is None:
            continue
        for j, item in zip(bidx, res):
            seqs[keys[j]] = item.get("seq", "")
        if (bi + 1) % 50 == 0 or bi == len(batches) - 1:
            print(f"      批 {bi+1}/{len(batches)}  已获取 {len(seqs):,}", flush=True)
        time.sleep(0.05)
    json.dump(seqs, open(CACHE, "w"))
    return seqs


print("[2] 从 Ensembl GRCh37 获取序列（+13bp 侧翼）...", flush=True)
seqs = run_fetch(ir)


# ---------------------------------------------------------------- 3. 特征
print("[3] 计算序列特征 ...", flush=True)
rows = []
drop = 0
for k, r in ir.iterrows():
    raw = seqs.get(k, "")
    if not raw:
        drop += 1
        continue
    ilen = int(r["ilen"])
    raw = raw.upper()
    if len(raw) != ilen + 6:          # 侧翼被截断跳过
        drop += 1
        continue
    if r["strand"] == "-":
        raw = revcomp(raw)
    intron = raw[3:3 + ilen]
    if len(intron) < 20:
        drop += 1
        continue
    gc = (intron.count("G") + intron.count("C")) / len(intron)
    ss5, ss3 = intron[:2], intron[-2:]
    tail = intron[-25:]
    pt20 = intron[-20:]
    y20 = (pt20.count("C") + pt20.count("T")) / 20
    maxy, cur = 0, 0
    for ch in tail:
        cur = cur + 1 if ch in "CT" else 0
        maxy = max(maxy, cur)
    # 5'ss 共识匹配：exon[-3..-1] = MAG，intron[+1..+6] = GTRAGT
    context5 = raw[:3] + intron[:6]          # 9-mer
    rows.append({
        "event": k, "gene": r["GeneSymbol"], "chrom": r["chrom"], "strand": r["strand"],
        "ilen": ilen, "e1len": r["e1len"], "e2len": r["e2len"],
        "gc": gc, "ss5": ss5, "ss3": ss3,
        "ss5_class": "GT" if ss5 == "GT" else ("GC" if ss5 == "GC" else ("AT" if ss5 == "AT" else "other")),
        "y20": y20, "maxY": maxy,
        "context5": context5,
        "MedianPSI": r["MedianPSI"], "stdPSI": r["stdPSI"], "pctSamples": r["pctSamples"],
        "novel": r.get("SpliceNovel", np.nan),
        "n_iso": len(str(r.get("SpliceInIsoName", "")).split(",")) if pd.notna(r.get("SpliceInIsoName")) else np.nan,
    })
F = pd.DataFrame(rows)
print(f"    可用于分析 {len(F):,}｜丢弃 {drop:,}（缺序列/过短/侧翼截断）")
F.to_csv(os.path.join(OUT, "G3_intron_features.tsv"), sep="\t", index=False)


# ---------------------------------------------------------------- 4. 分类
print("\n[4] IR 事件三类：Retained / Spliced / Variable")
print("    定义：Variable = stdPSI ≥ 0.10（高变异）；其余按 MedianPSI 三分位两端取 Spliced / Retained")
F = F.dropna(subset=["MedianPSI", "stdPSI"]).copy()
F["stdPSI"] = pd.to_numeric(F["stdPSI"], errors="coerce")
F["MedianPSI"] = pd.to_numeric(F["MedianPSI"], errors="coerce")
F = F.dropna(subset=["MedianPSI", "stdPSI"])
VAR_CUT = 0.10
grp = pd.Series("Intermediate", index=F.index)
grp[F["stdPSI"] >= VAR_CUT] = "Variable"
rest = F[F["stdPSI"] < VAR_CUT]
p = rest["MedianPSI"].quantile([1 / 3, 2 / 3]).values
grp[(F["stdPSI"] < VAR_CUT) & (F["MedianPSI"] <= p[0])] = "Spliced"
grp[(F["stdPSI"] < VAR_CUT) & (F["MedianPSI"] >= p[1])] = "Retained"
F["grp"] = grp
print(F["grp"].value_counts().to_string())
print(f"    MedianPSI 三分位切点：{p[0]:.3f} / {p[1]:.3f}")

from scipy import stats

NUM = ["ilen", "gc", "e1len", "e2len", "y20", "maxY", "n_iso"]
print("\n[5] 三类事件的特征对比（中位数 / 均值）")
print(f"{'feature':10s} {'Spliced':>14s} {'Retained':>14s} {'Variable':>14s} {'KW p':>12s} {'S vs R p':>12s}")
kw_all = {}
for c in NUM:
    g = [F.loc[F["grp"] == k, c].dropna() for k in ["Spliced", "Retained", "Variable"]]
    try:
        kw = stats.kruskal(*g)
    except Exception:
        kw = type("o", (), {"pvalue": np.nan})
    sr = stats.mannwhitneyu(g[0], g[1], alternative="two-sided")
    print(f"{c:10s} {g[0].median():14.3f} {g[1].median():14.3f} {g[2].median():14.3f} "
          f"{kw.pvalue:12.2e} {sr.pvalue:12.2e}")
    kw_all[c] = {"kw_p": kw.pvalue, "S_vs_R_p": sr.pvalue,
                 "med": [round(float(x.median()), 4) for x in g],
                 "mean": [round(float(x.mean()), 4) for x in g]}

print("\n[6] 剪接位点类别构成（%）")
ct = pd.crosstab(F["grp"], F["ss5_class"], normalize="index").mul(100).round(2)
print("  5' 剪接位点：")
print(ct.reindex(["Spliced", "Retained", "Variable"]).to_string())
ct2 = pd.crosstab(F["grp"], F["ss3"], normalize="index").mul(100).round(2)
print("\n  3' 剪接位点：")
print(ct2.reindex(["Spliced", "Retained", "Variable"])[["AG", "AC"]].to_string())
tab2 = stats.chi2_contingency(pd.crosstab(F.loc[F["grp"].isin(["Spliced", "Retained"]), "grp"],
                                          F.loc[F["grp"].isin(["Spliced", "Retained"]), "ss5_class"])[["GT", "GC", "AT", "other"]].values)
print(f"  非经典 5'ss 比例 Spliced vs Retained：卡方 p = {tab2[1]:.3e}")

# 长度匹配条件分析：长度是最重要的混杂
print("\n[7] 长度匹配后的条件分析（按内含子长度十分位分层）")
F["dlen"] = pd.qcut(F["ilen"], 10, labels=False, duplicates="drop")
res_cond = []
for c in ["gc", "y20", "maxY", "e1len", "e2len"]:
    S, R = F[F["grp"] == "Spliced"], F[F["grp"] == "Retained"]
    ds, dr, ws = [], [], []
    for d in sorted(set(F["dlen"].dropna())):
        a = S.loc[S["dlen"] == d, c].dropna(); b = R.loc[R["dlen"] == d, c].dropna()
        if len(a) < 30 or len(b) < 30:
            continue
        mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        w = len(a) * len(b) / (len(a) + len(b))
        ds.append(d); dr.append(mw.pvalue); ws.append(w)
        res_cond.append({"feature": c, "decile": int(d), "n_S": len(a), "n_R": len(b),
                         "med_S": round(a.median(), 4), "med_R": round(b.median(), 4),
                         "p": mw.pvalue, "w": w})
        if (d + 1) % 10 == 0 and c == "gc":
            pass
    # Stouffer 合并
    from scipy.stats import norm
    z = sum(norm.isf(p / 2) * (-1 if 0 else 1) * np.sqrt(w) for p, w in zip(dr, ws)) / np.sqrt(sum(ws))
    comb = 2 * norm.sf(abs(z))
    # 方向
    print(f"  {c:8s} 分层合并 Z = {z:+7.2f}  合并 p = {comb:.3e}  （层数 {len(ds)}）")
    for r0 in res_cond[-len(ds):]:
        pass
pd.DataFrame(res_cond).to_csv(os.path.join(OUT, "G3_length_matched.tsv"), sep="\t", index=False)

# 逻辑回归：预测 Retained(高PSI) vs Spliced
print("\n[8] 多因素逻辑回归：Retained (vs Spliced)")
try:
    import statsmodels.api as sm
    d = F[F["grp"].isin(["Spliced", "Retained"])].dropna(subset=NUM).copy()
    d["y"] = (d["grp"] == "Retained").astype(int)
    Xcols = ["ilen", "gc", "e1len", "e2len", "y20", "maxY", "n_iso"]
    Xraw = d[Xcols].fillna(d[Xcols].median())
    X = sm.add_constant(Xraw.copy())
    for c in ["ilen", "e1len", "e2len"]:
        X[c] = np.log10(X[c].clip(lower=1))
    fit = sm.Logit(d["y"], X).fit(disp=0)
    sd_model = X[Xcols].std()                      # 模型尺度的 SD
    sd_orig  = Xraw.std()                          # 原始尺度 SD（供注释）
    summ = pd.DataFrame({
        "coef": fit.params, "se": fit.bse, "z": fit.tvalues, "p": fit.pvalues,
        "OR_per_modelSD": np.exp(fit.params[Xcols] * sd_model),
    })
    summ.loc["const", "OR_per_modelSD"] = np.nan
    print(summ.round(4).to_string())
    summ.to_csv(os.path.join(OUT, "G3_logit.tsv"), sep="\t")
    print(f"  n = {len(d)}  pseudo-R² = {fit.prsquared:.4f}")
except Exception as e:
    print("  逻辑回归失败：", e)

json.dump({"kw": kw_all, "deciles": len(res_cond)}, open(os.path.join(OUT, "G3_summary.json"), "w"), indent=2)
print(f"\n输出：{OUT}/G3_intron_features.tsv  G3_length_matched.tsv  G3_logit.tsv")
print("完成")
