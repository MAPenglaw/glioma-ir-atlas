#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""19 — 年龄阳性对照 + burden 校正年龄后的多变量 Cox。

目的：证明 Kahles 2018 重跑版的「IR burden 非预后」阴性是真的（不是模型没效力），
用年龄作为阳性对照：年龄在 GBM/LGG 都应显著预后，而 burden 不显著。
另跑 burden + 年龄 双变量，看校正后 burden 是否仍不显著。
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import stats

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
OUT = os.path.join(ROOT, "ri_atlas", "output_km")
CLIN = os.path.join(ROOT, "data", "TCGA_GBMLGG_clinical.tsv")


def barcode(c):
    return str(c).strip().strip('"')[:15]


def cox_ph(time, event, *xs):
    """单变量或多变量 Cox（xs 为标准化后的列）。返回 (hr, hr_lo, hr_hi, p)。"""
    X = np.column_stack(xs)
    n = X.shape[0]
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    order = np.argsort(time)
    time = time[order]; event = event[order]; X = X[order]

    def neglog(beta):
        b = np.asarray(beta, dtype=float)
        expb = np.exp(X @ b)
        cum = np.cumsum(expb[::-1])[::-1]
        pl = np.sum(event * (X @ b - np.log(cum + 1e-300)))
        return -pl

    k = X.shape[1]
    res = minimize(neglog, np.zeros(k), method="Nelder-Mead")
    beta = res.x
    # 数值 Hessian
    eps = 1e-4
    hess = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            b1 = beta.copy(); b1[i] += eps; b1[j] += eps
            b2 = beta.copy(); b2[i] += eps; b2[j] -= eps
            b3 = beta.copy(); b3[i] -= eps; b3[j] += eps
            b4 = beta.copy(); b4[i] -= eps; b4[j] -= eps
            hess[i, j] = (neglog(b1) - neglog(b2) - neglog(b3) + neglog(b4)) / (4 * eps ** 2)
    cov = np.linalg.inv(hess + 1e-6 * np.eye(k))
    out = []
    for i in range(k):
        se = np.sqrt(max(cov[i, i], 1e-9))
        z = beta[i] / se
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        hr = np.exp(beta[i])
        hr_lo = np.exp(beta[i] - 1.96 * se)
        hr_hi = np.exp(beta[i] + 1.96 * se)
        out.append((hr, hr_lo, hr_hi, p))
    return out


def run(cancer):
    outdir = os.path.join(OUT, cancer)
    burden = pd.read_csv(os.path.join(outdir, "burden.tsv"), sep="\t", index_col=0)["RI_burden"]
    burden.index = [barcode(i) for i in burden.index]
    burden = burden.groupby(level=0).mean()

    surv = pd.read_csv(os.path.join(ROOT, "ri_atlas", "data", f"{cancer}_survival.txt"), sep="\t")
    surv["sample"] = [barcode(s) for s in surv["sample"]]
    surv = surv.drop_duplicates(subset="sample").set_index("sample")

    clin = pd.read_csv(CLIN, sep="\t", usecols=["sampleID", "age_at_initial_pathologic_diagnosis"])
    clin["sample"] = [barcode(s) for s in clin["sampleID"]]
    clin = clin.drop_duplicates(subset="sample").set_index("sample")
    age = pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce")

    common = burden.index.intersection(surv.index).intersection(age.dropna().index)
    b = burden.reindex(common).astype(float)
    a = age.reindex(common).astype(float)
    os_time = surv.reindex(common)["OS.time"].astype(float)
    os_event = surv.reindex(common)["OS"].astype(float)
    print(f"[{cancer}] 交集 n = {len(common)}（burden ∩ 生存 ∩ 年龄）")

    # 年龄单变量（阳性对照），年龄按每10岁，HR 更直观
    a10 = a / 10.0
    r = cox_ph(os_time, os_event, (a10 - a10.mean()) / (a10.std() + 1e-12))
    hr, lo, hi, p = r[0]
    print(f"  年龄(每10岁) vs OS: HR={hr:.3f} [{lo:.3f}-{hi:.3f}] p={p:.3g}  ← 阳性对照")

    # burden 单变量（每 SD）
    bz = (b - b.mean()) / (b.std() + 1e-12)
    r = cox_ph(os_time, os_event, bz)
    hr, lo, hi, p = r[0]
    print(f"  burden(每SD) vs OS: HR={hr:.3f} [{lo:.3f}-{hi:.3f}] p={p:.3g}")

    # burden + 年龄 双变量
    az = (a10 - a10.mean()) / (a10.std() + 1e-12)
    r = cox_ph(os_time, os_event, bz, az)
    print(f"  burden+年龄: burden HR={r[0][0]:.3f} p={r[0][3]:.3g} | 年龄 HR={r[1][0]:.3f} p={r[1][3]:.3g}")

    return common, b, a, os_time, os_event


if __name__ == "__main__":
    for c in ("GBM", "LGG"):
        run(c)
        print()
