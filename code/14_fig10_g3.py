#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G3 出图：Fig10_G3_intron_features（PDF 矢量 + PNG）"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/output"
FIG = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/figures"

C = {"Spliced": "#6B8F71", "Retained": "#C1553B", "Variable": "#3D7EA6"}
order = ["Spliced", "Retained", "Variable"]
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#8A8378",
                     "axes.linewidth": 0.9, "axes.labelsize": 9.5,
                     "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "axes.titlesize": 10.5, "axes.titleweight": "bold"})
BG = "#FCFAF6"

F = pd.read_csv(f"{OUT}/G3_intron_features.tsv", sep="\t")
F["stdPSI"] = pd.to_numeric(F["stdPSI"], errors="coerce")
F["MedianPSI"] = pd.to_numeric(F["MedianPSI"], errors="coerce")
F = F.dropna(subset=["stdPSI", "MedianPSI"])
grp = pd.Series("Intermediate", index=F.index)
grp[F["stdPSI"] >= 0.10] = "Variable"
r = F[F["stdPSI"] < 0.10]
q = r["MedianPSI"].quantile([1/3, 2/3]).values
grp[(F["stdPSI"] < 0.10) & (F["MedianPSI"] <= q[0])] = "Spliced"
grp[(F["stdPSI"] < 0.10) & (F["MedianPSI"] >= q[1])] = "Retained"
F["grp"] = grp
lm = pd.read_csv(f"{OUT}/G3_length_matched.tsv", sep="\t")

fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6), facecolor=BG)
fig.subplots_adjust(hspace=0.42, wspace=0.44, left=0.06, right=0.97, top=0.90, bottom=0.08)

def style(ax, title):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(title, loc="left", pad=8)
    ax.grid(axis="y", color="#E4DED2", lw=0.7, zorder=0)

# ---------------- A 内含子长度 ----------------
ax = axes[0, 0]
data = [np.log10(F.loc[F["grp"] == g, "ilen"].clip(lower=1)) for g in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55, showfliers=False,
                medianprops=dict(color="white", lw=1.6))
for p, g in zip(bp["boxes"], order):
    p.set_facecolor(C[g]); p.set_alpha(0.85); p.set_edgecolor("none")
ax.set_ylabel("Intron length (bp, log10)")
ax.set_yticks([2, 2.5, 3, 3.5]); ax.set_yticklabels(["100", "316", "1k", "3.2k"])
style(ax, "A  Retained introns are SHORTER")
ax.text(0.03, 0.04, "median 264 / 194 / 374 bp", transform=ax.transAxes, fontsize=8.2,
        color="#5A544B", style="italic")

# ---------------- B 5'ss 构成 ----------------
ax = axes[0, 1]
ct = pd.crosstab(F["grp"], F["ss5_class"], normalize="index").mul(100)
ct = ct.reindex(order)[["GT", "GC", "AT", "other"]]
bottom = np.zeros(3)
cols5 = {"GT": "#8AA98C", "GC": "#C1553B", "AT": "#7C5CB0", "other": "#B8AFA2"}
for k in ct.columns:
    ax.bar(range(3), ct[k], bottom=bottom, color=cols5[k], width=0.62, label=k)
    bottom += ct[k].values
for i, g in enumerate(order):
    nongt = 100 - ct.loc[g, "GT"]
    if nongt > 0.5:
        ax.text(i, 101, f"{nongt:.1f}%", ha="center", fontsize=8.6, color=C["Retained"],
                weight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels(order)
ax.set_ylim(0, 132); ax.set_ylabel("% of events")
ax.legend(ncol=4, frameon=False, fontsize=8, loc="lower center",
          bbox_to_anchor=(0.5, -0.28), title="5' splice site", title_fontsize=8)
style(ax, "B  Non-canonical 5'ss enriched 20× in Retained")
ax.text(0.5, 0.965, "GC-AG: 17.2% (Retained) vs 0.78% (Spliced) · χ² p = 3.9e-206",
        transform=ax.transAxes, ha="center", fontsize=8.3, color="#C1553B", weight="bold")

# ---------------- C poly-Y ----------------
ax = axes[0, 2]
data = [F.loc[F["grp"] == g, "y20"].dropna() for g in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55, showfliers=False,
                medianprops=dict(color="white", lw=1.6))
for p, g in zip(bp["boxes"], order):
    p.set_facecolor(C[g]); p.set_alpha(0.85); p.set_edgecolor("none")
ax.set_ylabel("Poly-Y content, last 20 nt")
style(ax, "C  Weak polypyrimidine tract in Retained")
ax.text(0.03, 0.92, "length-matched merged p < 1e-300", transform=ax.transAxes, fontsize=8.2,
        color="#5A544B", style="italic")

# ---------------- D GC ----------------
ax = axes[1, 0]
data = [F.loc[F["grp"] == g, "gc"].dropna() for g in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55, showfliers=False,
                medianprops=dict(color="white", lw=1.6))
for p, g in zip(bp["boxes"], order):
    p.set_facecolor(C[g]); p.set_alpha(0.85); p.set_edgecolor("none")
ax.set_ylabel("Intron GC content")
ax.set_ylim(0.30, 0.75)
style(ax, "D  Higher GC in Retained (univariate)")
ax.text(0.03, 0.92, "not independent in multivariable model", transform=ax.transAxes,
        fontsize=8.2, color="#5A544B", style="italic")

# ---------------- E 长度分层条件中位数 ----------------
ax = axes[1, 1]
for feat, col, lab in [("y20", "#7C5CB0", "poly-Y (y20)"),
                       ("gc", "#C1553B", "GC content")]:
    sub = lm[lm["feature"] == feat].sort_values("decile")
    ax.plot(sub["decile"] + 1, sub["med_R"] - sub["med_S"], "o-", color=col, lw=1.6,
            ms=4.5, label=f"{lab} (R − S)")
ax.axhline(0, color="#8A8378", lw=0.9, ls="--")
ax.set_ylim(-0.32, 0.115)
ax.set_xlabel("Intron-length decile (1 = shortest)")
ax.set_ylabel("Median difference (Retained − Spliced)")
ax.legend(frameon=False, fontsize=8.0, loc="lower right", bbox_to_anchor=(0.99, 0.04))
style(ax, "E  Differences persist at every length stratum")
ax.text(0.5, 0.94, "Stouffer merged p:  y20 < 1e-300 · GC = 2.2e-14", transform=ax.transAxes,
        fontsize=8.2, color="#5A544B", style="italic", ha="center", va="top")

# ---------------- F 逻辑回归森林 ----------------
ax = axes[1, 2]
lg = pd.read_csv(f"{OUT}/G3_logit.tsv", sep="\t", index_col=0).drop(index="const")
lg = lg.drop(index="maxY")          # 与 y20 共线，主模型不含（敏感性见报告）
lg = lg.loc[["e1len", "e2len", "n_iso", "gc", "ilen", "y20"]]
lab = {"e1len": "Upstream exon len (log10)", "e2len": "Downstream exon len (log10)",
       "n_iso": "Isoform count", "gc": "GC content", "ilen": "Intron len (log10)",
       "y20": "Poly-Y content"}
y = np.arange(len(lg))[::-1]
orv = lg["OR_per_modelSD"].values
se = lg["se"].values
lo, hi = np.exp(np.log(orv) - 1.96 * se / orv * orv), np.exp(np.log(orv) + 1.96 * se / orv * orv)
ax.axvline(1, color="#8A8378", lw=0.9, ls="--")
for i, (g, o) in enumerate(zip(lg.index, orv)):
    col = "#C1553B" if o > 1 else "#3D7EA6"
    yi = y[i]
    ax.plot([lo[i], hi[i]], [yi, yi], color=col, lw=1.8)
    ax.plot(o, yi, "o", color=col, ms=6)
ax.set_yticks(y); ax.set_yticklabels([lab[g] for g in lg.index], fontsize=8.6)
ax.set_xscale("log"); ax.set_xlim(0.05, 30)
ax.set_xlabel("OR per model-SD (Retained vs Spliced)")
style(ax, "F  Multivariable logistic regression")
ax.text(0.03, 0.10, "n = 2,832 · pseudo-R² = 0.37\n(maxY omitted: collinear with y20)",
        transform=ax.transAxes, fontsize=8.2, color="#5A544B", style="italic")

fig.suptitle("Genomic features of intron-retention events in TCGA-GBM (n = 20,263 IR events)",
             fontsize=12.5, weight="bold", y=0.975, color="#2E2B28")
for ext in ("pdf", "png"):
    p = os.path.join(FIG, f"Fig10_G3_intron_features.{ext}")
    fig.savefig(p, dpi=300 if ext == "png" else None, facecolor=BG,
                bbox_inches="tight", pad_inches=0.25)
    print("写出:", p)
plt.close(fig)
