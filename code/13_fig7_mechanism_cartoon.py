#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G5：机制卡通模型图（Fig 7）
输出真矢量 PDF + 300dpi PNG，全英文标注，可直接进论文。
所有数字均为本机实测（见 G1G2_执行报告.md / RI图谱_方案B执行报告.md）。
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.patches import Rectangle
import numpy as np

OUT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas/figures"
os.makedirs(OUT, exist_ok=True)

# ---------------- palette ----------------
C_DOWN  = "#C1553B"      # 下降 / 不足
C_UP    = "#3D7EA6"      # 上升 / 代偿
C_NEUT  = "#88796A"
C_TEXT  = "#2E2B28"
C_BG    = "#FCFAF6"
C_PANEL = "#F4EFE6"
C_LINE  = "#B8AFA2"
C_OK    = "#6B8F71"      # 正常态绿
grN     = dict(fc=C_PANEL, ec=C_LINE, lw=1.1)

fig = plt.figure(figsize=(13.6, 7.4), facecolor=C_BG)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, **kw):
    d = dict(fc=C_PANEL, ec=C_LINE, lw=1.1, boxstyle="round,pad=0.35,rounding_size=1.1")
    d.update(kw)
    ax.add_patch(FancyBboxPatch((x, y), w, h, **d))

def txt(x, y, s, size=9.5, col=C_TEXT, weight="normal", ha="center", va="center",
        style="normal", family="DejaVu Sans"):
    ax.text(x, y, s, fontsize=size, color=col, weight=weight, ha=ha, va=va,
            style=style, family=family, zorder=20)

def arrow(p1, p2, col=C_NEUT, lw=1.8, ls="-", rad=0.0, head=7, arrowstyle="-|>"):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=arrowstyle, mutation_scale=head,
                                 color=col, lw=lw, linestyle=ls, zorder=10,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=1.5, shrinkB=1.5))

def header(x, y, tag, title):
    ax.text(x, y, tag, fontsize=15, color=C_DOWN, weight="bold", family="DejaVu Sans",
            ha="left", va="center", zorder=20)
    ax.text(x + 3.2, y, title, fontsize=12.5, color=C_TEXT, weight="bold",
            family="DejaVu Sans", ha="left", va="center", zorder=20)

# =====================================================================
# 顶部标题条
# =====================================================================
box(3, 92.5, 94, 6.2, fc="#EFE7DA", ec=C_LINE, lw=1.0,
    boxstyle="round,pad=0.2,rounding_size=0.9")
txt(50, 95.6, "Proposed model: spliceosome assembly-core insufficiency drives a global intron-retention burden in glioma",
    size=11.5, weight="bold")

# =====================================================================
# A. Homeostasis  (x 4 → 40)
# =====================================================================
header(5.5, 87.5, "A", "Spliceosome homeostasis")

# 基因表达层面
box(5.5, 75.5, 34, 8.2, **grN)
txt(22.5, 81.4, "Sm-core / Sm-ring subunits  (SNRPB/D/E/F/G …)", size=9.2, weight="bold")
txt(22.5, 78.3, "abundant  →  U1 / U2 / U4 / U5 snRNP pool sufficient", size=8.6, col=C_NEUT)
arrow((22.5, 75.3), (22.5, 71.8), C_OK, lw=2.0, head=9)

box(5.5, 62.5, 34, 9.0, fc="#E7EDE6", ec="#A9BCA6", lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(22.5, 68.9, "Spliceosome assembles efficiently", size=9.4, weight="bold", col="#4F6B53")
txt(22.5, 65.6, "B / Bact / C complexes form normally", size=8.6, col=C_NEUT)
arrow((22.5, 62.2), (22.5, 58.6), C_OK, lw=2.0, head=9)

# pre-mRNA → mature
box(5.5, 49.5, 34, 8.6, fc="#EFE7DA", ec=C_LINE, lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(22.5, 55.6, "Pre-mRNA", size=9.2, weight="bold")
# 内含子被剪掉示意
ax.add_patch(Rectangle((10.5, 50.9), 24, 1.5, fc="#6B8F71", ec="none", zorder=15))
txt(22.5, 51.6, "intron excised", size=7.8, col="white", weight="bold")
arrow((22.5, 49.2), (22.5, 45.8), C_OK, lw=2.0, head=9)

box(5.5, 37.5, 34, 8.0, fc="#E7EDE6", ec="#A9BCA6", lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(22.5, 43.4, "Mature mRNA  →  protein output", size=9.2, weight="bold", col="#4F6B53")
txt(22.5, 40.1, "low global intron-retention burden", size=8.6, col=C_NEUT)

txt(22.5, 33.6, "Intron retention is a minor,\nstochastic event", size=8.8, col=C_OK, weight="bold")

# =====================================================================
# B. Insufficiency  (x 60 → 96)
# =====================================================================
header(60.5, 87.5, "B", "Assembly-core insufficiency  (observed in glioma)")

box(60.5, 75.5, 34, 8.2, fc="#F7ECE7", ec="#CFB2A3", lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(77.5, 81.4, "Sm-core subunits DOWN", size=9.4, weight="bold", col=C_DOWN)
txt(77.5, 78.3, "7/7 members negative · median ρ = −0.294 · perm p = 0.022", size=8.4, col=C_DOWN)
arrow((77.5, 75.3), (77.5, 71.8), C_DOWN, lw=2.0, head=9)

box(60.5, 62.5, 34, 9.0, fc="#F7ECE7", ec="#CFB2A3", lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(77.5, 68.9, "snRNP supply bottleneck", size=9.4, weight="bold", col=C_DOWN)
txt(77.5, 65.6, "NineTeen / Bact 10/12 negative · p = 0.059", size=8.4, col=C_DOWN)
arrow((77.5, 62.2), (77.5, 58.6), C_DOWN, lw=2.0, head=9)

# retained
box(60.5, 49.5, 34, 8.6, fc="#EFE7DA", ec=C_LINE, lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(77.5, 55.6, "Pre-mRNA with RETAINED intron", size=9.2, weight="bold", col=C_DOWN)
ax.add_patch(Rectangle((65.0, 50.9), 25, 1.5, fc=C_DOWN, ec="none", zorder=15))
txt(77.5, 51.6, "intron retained", size=7.8, col="white", weight="bold")
arrow((77.5, 49.2), (77.5, 45.8), C_DOWN, lw=2.0, head=9)

box(60.5, 37.5, 34, 8.0, fc="#F7ECE7", ec="#CFB2A3", lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
txt(77.5, 43.4, "NMD / nuclear decay  →  mRNA ↓", size=9.2, weight="bold", col=C_DOWN)
txt(77.5, 40.1, "20,365 IR events (29.6% of all events)", size=8.6, col=C_DOWN)

# 代偿支
arrow((64.0, 45.4), (52.0, 45.4), C_UP, lw=1.7, head=8)
box(41.0, 41.9, 11.5, 7.2, fc="#E7EEF5", ec="#A9BED1", lw=1.1,
    boxstyle="round,pad=0.3,rounding_size=0.9")
txt(46.75, 47.6, "Regulatory layer", size=8.3, weight="bold", col=C_UP)
txt(46.75, 45.2, "SR proteins UP", size=8.4, weight="bold", col=C_UP)
txt(46.75, 42.6, "ρ = +0.361", size=8.2, col=C_UP)
arrow((46.75, 41.6), (46.75, 38.6), C_UP, lw=1.5, head=7)
txt(46.75, 36.4, "compensatory,\ninsufficient", size=7.8, col=C_UP, style="italic")

# SNRPC 标签：贴在 B 侧
box(60.5, 26.4, 34, 6.0, fc="#F2E3D9", ec="#C98F64", lw=1.2,
    boxstyle="round,pad=0.3,rounding_size=0.9")
txt(77.5, 30.6, "SNRPC (U1-C)", size=8.8, weight="bold", col=C_DOWN)
txt(77.5, 27.6, "ρ = −0.320 · BH q = 2.1e-04", size=8.3, col=C_DOWN)

txt(22.5, 26.6, "External causal anchor (GSE227187):\nSNRPC knockdown → 1,400 IR events → 21/30 (70%) with\nmRNA down-regulation · OXPHOS 77% down (p = 5.8e-17)",
    size=8.3, col=C_TEXT, weight="normal")
arrow((30.0, 24.4), (39.0, 29.6), C_NEUT, lw=1.3, head=7, ls=(0, (4, 3)))
txt(44.5, 25.0, "SUPPORTIVE, not\nglioma-intrinsic\n(red line R7)", size=7.6, col="#9A9086", style="italic")

# =====================================================================
# C. 连续轴 + 结论
# =====================================================================
box(3, 1.5, 94, 20.5, fc="#F7F2E9", ec=C_LINE, lw=1.1,
    boxstyle="round,pad=0.35,rounding_size=1.1")
header(5.5, 19.6, "C", "It is a continuous burden axis — and it is NOT a prognostic biomarker")

# 梯度条
g_ax = fig.add_axes([0.075, 0.128, 0.60, 0.032])
grad = np.linspace(0, 1, 256).reshape(1, -1)
g_ax.imshow(grad, aspect="auto", cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
    "gb", ["#EDE6D8", "#D8C0A6", "#C98F64", "#C1553B"]))
g_ax.set_xticks([]); g_ax.set_yticks([])
for s in g_ax.spines.values():
    s.set_color(C_LINE); s.set_linewidth(1.0)
g_ax.text(-0.018, 0.5, "low", transform=g_ax.transAxes, fontsize=8.6, color=C_NEUT,
          ha="right", va="center")
g_ax.text(1.018, 0.5, "high", transform=g_ax.transAxes, fontsize=8.6, color=C_DOWN,
          ha="left", va="center")
txt(33.5, 9.9, "Global intron-retention burden (low → high) · continuous · PC1 = 25.3% variance · no discrete subtypes",
    size=8.5, col=C_TEXT)

# 右侧结论块
box(76.0, 6.4, 19.5, 13.0, fc="#F7ECE7", ec="#CFB2A3", lw=1.1,
    boxstyle="round,pad=0.3,rounding_size=0.9")
txt(85.75, 17.4, "NOT prognostic", size=10.2, weight="bold", col=C_DOWN)
txt(85.75, 14.6, "OS HR/SD = 0.976 (p = 0.78)", size=8.4, col=C_TEXT)
txt(85.75, 12.5, "0 / 16 metrics × 2 endpoints", size=8.4, col=C_TEXT)
txt(85.75, 10.4, "CGGA after adjustment: p = 0.54", size=8.4, col=C_TEXT)
txt(85.75, 8.3, "positive control age: p = 6e-4", size=8.2, col=C_NEUT, style="italic")

txt(5.5, 6.2, "Localisation:  burden-axis genes peak in malignant cells",
    size=8.8, col=C_TEXT, ha="left", weight="bold")
txt(5.5, 3.9, "GSE182109 · 264,523 cells · 0.244 vs 0.189 · p < 1e-300",
    size=8.3, col=C_TEXT, ha="left")
txt(43.5, 6.2, "Decoupled from confounders",
    size=8.8, col=C_TEXT, ha="left", weight="bold")
txt(43.5, 3.9, "proliferation ρ = −0.045 (p = 0.59) · purity ρ = −0.014 · stroma ρ = −0.011",
    size=8.3, col=C_TEXT, ha="left")

# 中间竖向分隔线
ax.plot([50, 50], [28.0, 90.0], color=C_LINE, lw=0.9, ls=(0, (5, 4)), zorder=5)

for ext in ("pdf", "png"):
    p = os.path.join(OUT, f"Fig7_mechanism_model.{ext}")
    fig.savefig(p, dpi=300 if ext == "png" else None, facecolor=C_BG,
                bbox_inches="tight", pad_inches=0.25)
    print("写出:", p)
plt.close(fig)
