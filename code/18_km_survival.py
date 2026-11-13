#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""18 — Kahles 2018 数据的生存分析：burden 是否预后（核心阴性结论验证）。

用法: python 18_km_survival.py GBM | LGG
单变量 Cox（手动 partial likelihood）+ log-rank（burden 中位数分组）。
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
OUT = os.path.join(ROOT, "ri_atlas", "output_km")


def barcode(c):
    return str(c).strip().strip('"')[:15]


def cox_ph(time, event, x):
    """单变量 Cox 比例风险。返回 (beta, se, z, p, hr, hr_lo, hr_hi)。"""
    # x 标准化
    x = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
    order = np.argsort(time)
    time = time[order]; event = event[order]; x = x[order]
    # 负 log partial likelihood
    def neglog(beta):
        b = beta[0]
        expb = np.exp(b * x)
        # 风险集累计和（从后往前）
        cum = np.cumsum(expb[::-1])[::-1]
        # 只对 event=1 的样本贡献
        pl = np.sum(event * (b * x - np.log(cum + 1e-300)))
        return -pl
    res = minimize(neglog, [0.0], method="Nelder-Mead")
    beta = res.x[0]
    # Hessian via 数值二阶导
    eps = 1e-4
    d2 = (neglog([beta + eps]) - 2 * neglog([beta]) + neglog([beta - eps])) / eps**2
    se = np.sqrt(1.0 / max(d2, 1e-6))
    z = beta / se
    p = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z)))
    hr = np.exp(beta)
    hr_lo = np.exp(beta - 1.96 * se)
    hr_hi = np.exp(beta + 1.96 * se)
    return beta, se, z, p, hr, hr_lo, hr_hi


def logrank(time, event, group):
    """两分组 log-rank。返回 (chi2, p)。"""
    import numpy as np
    u = np.unique(time)
    u = u[u >= 0]
    O1 = E1 = V = 0.0
    for t in u:
        at = time == t
        d = event[at].sum()
        n = at.sum()
        if n == 0:
            continue
        n1 = (group[at] == 1).sum()
        d1 = ((event[at] > 0.5) & (group[at] == 1)).sum()
        if n > 1 and d > 0:
            E1 += d * n1 / n
            V += d * (n - d) * n1 * (n - n1) / (n**2 * (n - 1))
        O1 += d1
    chi2 = (O1 - E1) ** 2 / (V + 1e-9)
    p = 1 - __import__("scipy").stats.chi2.cdf(chi2, 1)
    return chi2, p


def main():
    cancer = sys.argv[1] if len(sys.argv) > 1 else "GBM"
    outdir = os.path.join(OUT, cancer)
    burden = pd.read_csv(os.path.join(outdir, "burden.tsv"), sep="\t", index_col=0)["RI_burden"]
    burden.index = [barcode(i) for i in burden.index]
    if burden.index.duplicated().any():
        burden = burden.groupby(level=0).mean()

    surv = pd.read_csv(os.path.join(ROOT, "ri_atlas", "data", f"{cancer}_survival.txt"), sep="\t")
    surv["sample"] = [barcode(s) for s in surv["sample"]]
    surv = surv.drop_duplicates(subset="sample").set_index("sample")

    common = burden.index.intersection(surv.index)
    print(f"[1] {cancer}: burden {len(burden)}，生存数据交集 {len(common)}")

    b = burden.reindex(common)
    rows = []
    for endpoint in ["OS", "DSS", "PFI"]:
        time = pd.to_numeric(surv.loc[common, f"{endpoint}.time"], errors="coerce").values
        event = pd.to_numeric(surv.loc[common, endpoint], errors="coerce").values
        ok = ~np.isnan(time) & ~np.isnan(event) & (time > 0)
        t, e, x = time[ok], event[ok], b.values[ok]
        if len(t) < 30:
            print(f"    {endpoint}: 有效样本不足 {len(t)}，跳过")
            continue
        beta, se, z, p, hr, lo, hi = cox_ph(t, e, x)
        med = np.median(x)
        grp = (x > med).astype(int)
        chi2, plr = logrank(t, e, grp)
        rows.append({"endpoint": endpoint, "n": int(ok.sum()), "HR": hr, "CI_lo": lo,
                     "CI_hi": hi, "p_cox": p, "p_logrank": plr})
        print(f"    {endpoint:<4} n={ok.sum():>4}  HR={hr:.3f} [{lo:.3f}–{hi:.3f}] "
              f"p_cox={p:.3f}  p_logrank={plr:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "survival.tsv"), sep="\t", index=False)
    print("[2] 完成 ->", os.path.join(outdir, "survival.tsv"))


if __name__ == "__main__":
    main()
