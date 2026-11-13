#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G2b — 剪接体亚复合体层面：哪一部分机器的表达与 RI 负荷相关？

背景：G2 发现 U1 panel 内部方向不一致（Sm 核心全负、SNRNP70 强正），
      说明「U1 snRNP 整体下降」的表述不准确，需要拆到亚复合体层面。
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = f"{BASE}/output"
N_PERM = 20000
SEED = 20260917

GROUPS = {
    "Sm 核心（U1/U2/U4/U5 共用）": ["SNRPB", "SNRPD1", "SNRPD2", "SNRPD3",
                                    "SNRPE", "SNRPF", "SNRPG"],
    "U1 特异亚基": ["SNRPC", "SNRPA", "SNRNP70"],
    "U2 / 3' 剪接位点识别": ["SF3A1", "SF3B1", "SF3B2", "SF3B3", "U2AF1", "U2AF2",
                             "RBM39", "PHF5A", "SF1"],
    "U4/U6.U5 tri-snRNP": ["SNRNP200", "SNRNP40", "PRPF8", "EFTUD2", "TXNL4A",
                            "PPIH", "PRPF6", "PRPF31", "PRPF3", "PRPF4",
                            "SNRNP27", "USP39", "SART1", "SNRNP48"],
    "NineTeen / Bact 复合体": ["ISY1", "ZMAT2", "CDC5L", "BCAS2", "BUD31",
                                "CRNKL1", "AQR", "XAB2", "PPIE", "SYF2",
                                "CWC15", "CTNNBL1", "SNW1", "PLRG1", "PRPF19"],
    "SR 蛋白家族": ["SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6",
                    "SRSF7", "SRSF8", "SRSF9", "SRSF10", "SRSF11", "SRSF12"],
    "hnRNP 家族": ["HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPF", "HNRNPH1",
                   "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPU", "PTBP1", "PTBP2"],
    "3' 端加工 / 输出": ["CPSF1", "CPSF2", "CPSF3", "CPSF4", "CPSF6", "CPSF7",
                         "CSTF1", "CSTF2", "CSTF3", "NUDT21", "PABPN1",
                         "ALYREF", "NXF1", "MAGOH", "CASC3"],
}


def main():
    df = pd.read_csv(f"{OUT}/G2_sf_vs_burden.tsv", sep="\t")
    d = dict(zip(df["gene"], df["rho"]))
    all_rho = df["rho"].dropna().values
    print(f"载入 {len(df)} 个剪接因子的结果")
    print("=" * 84)
    print(f"{'亚复合体':26s} {'n':>4s} {'负相关':>8s} {'中位ρ':>9s} "
          f"{'二项p':>10s} {'置换p':>9s}  成员（ρ）")
    print("=" * 84)

    rng = np.random.default_rng(SEED)
    rows = []
    for gname, gl in GROUPS.items():
        gs = [g for g in gl if g in d]
        if len(gs) < 3:
            continue
        rs = np.array([d[g] for g in gs])
        n_neg = int((rs < 0).sum())
        med = float(np.median(rs))
        bt = stats.binomtest(n_neg, len(rs), 0.5, alternative="two-sided")
        null = np.array([np.mean(rng.choice(all_rho, len(rs), replace=False) < 0)
                         for _ in range(N_PERM)])
        p_perm = (np.sum(np.abs(null - 0.5) >= abs(n_neg / len(rs) - 0.5)) + 1) / (N_PERM + 1)
        detail = ", ".join(f"{g}{d[g]:+.2f}" for g in gs)
        print(f"{gname:26s} {len(gs):4d} {n_neg/len(gs):7.1%} {med:+9.4f} "
              f"{bt.pvalue:10.3g} {p_perm:9.4g}  {detail}")
        rows.append({"group": gname, "n": len(gs), "n_neg": n_neg,
                     "frac_neg": n_neg / len(gs), "median_rho": med,
                     "p_binom": bt.pvalue, "p_perm": p_perm,
                     "mean_rho": float(np.mean(rs)), "members": "|".join(gs)})

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/G2b_complex_groups.tsv", sep="\t", index=False)
    print(f"\n已输出 {OUT}/G2b_complex_groups.tsv")

    # 显著性成员汇总
    sig = df[df["q"] < 0.05].sort_values("rho")
    print(f"\n全池 BH q<0.05 的基因：{len(sig)} 个（负 {int((sig['rho']<0).sum())}"
          f" / 正 {int((sig['rho']>0).sum())}）")


if __name__ == "__main__":
    main()
