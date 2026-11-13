#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""56 — 事件数匹配的下采样对照：验证 SpliceSeq vs SplAdder 的聚类指标差异
是否只是「高变异事件数量不同」的技术效应。

背景
----
55 号脚本用第二套管道 TCGASpliceSeq 复现，结果分化为两类：
    T1 负担轴一致性        ρ=0.960 (GBM) / 0.958 (LGG)   —— 完美复现
    T2 聚类诊断            silhouette 0.45/0.39 vs 0.148/0.124 —— 显著更高
                           LGG 的 GMM-BIC 甚至从 k1 翻转到 k2
但两套管道的事件规模差异巨大：
    SplAdder   高变异事件 GBM 1509 / LGG 2216（自 15 万+ IR 事件筛出）
    SpliceSeq  高变异事件 GBM  472 / LGG  300（自 ~2.3k/1.3k RI 事件筛出）
silhouette 对「样本数 / 特征数」比极度敏感：特征越少，同一次 KMeans 越容易
在低维投影里切出看似分离的簇。若不排除这个解释就写"第二条被推翻"，是错的。

做法
----
在 **SplAdder 自己的高变异事件池**里，随机下采样到 SpliceSeq 的事件数量
（GBM 472 / LGG 300），重复 20 次，看 silhouette / BIC 偏好 / PAC 的分布。
    - 若下采样后的 SplAdder 指标逼近 SpliceSeq 的值  → 差异是事件数效应，
      "无离散亚型"的结论未被推翻，反而需要如实披露这个敏感性。
    - 若下采样后仍远低于 SpliceSeq  → 差异是管道本身的真实差异，必须改写结论。

输出 output_secondpipeline/downsample_control.csv
"""
import os
import time

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

ROOT_B = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
ATLAS = os.path.join(ROOT_B, "ri_atlas")
OUT = os.path.join(ATLAS, "output_secondpipeline")
os.makedirs(OUT, exist_ok=True)

MIN_DETECT = 0.90
MIN_STD = 0.10
SEED = 20260917
N_REP = 10
N_ITER_PAC = 60
SUB_SAMPLE = 0.80
SUB_FEATURE = 0.80

# SpliceSeq 侧实测的高变异事件数（55 号脚本输出），就是要匹配的目标
TARGET = {"GBM": 472, "LGG": 300}
# SpliceSeq 侧实测值，作对照基准
SPLICESEQ_VAL = {
    "GBM": {"sil": 0.450, "pac": 0.521, "bic": "k1"},
    "LGG": {"sil": 0.394, "pac": 0.487, "bic": "k2"},
}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_spladder_hivar(cancer):
    psi = pd.read_csv(f"{ATLAS}/data/KM_IR_{cancer}.psi.gz", sep="\t",
                      index_col=0, compression="gzip")
    vals = psi.iloc[:, 4:].apply(pd.to_numeric, errors="coerce")
    vals.columns = [str(c).split(".")[0] for c in vals.columns]
    vals = vals.loc[:, ~vals.columns.duplicated()]
    detect = vals.notna().mean(axis=1)
    sd = vals.std(axis=1, skipna=True)
    hv = vals[(detect >= MIN_DETECT) & (sd >= MIN_STD)]
    X = hv.T
    X = X.apply(lambda c: c.fillna(c.mean()), axis=0)
    return X.dropna(axis=1, how="all")


def pac_fixed_old(Z, seed):
    n, p = Z.shape
    n_s = max(2, int(round(n * SUB_SAMPLE)))
    n_f = max(2, int(round(p * SUB_FEATURE)))
    rng = np.random.default_rng(seed)
    cnt = {}
    for it in range(N_ITER_PAC):
        si = rng.choice(n, n_s, replace=False)
        fi = rng.choice(p, n_f, replace=False)
        Zi = Z[np.ix_(si, fi)]
        lab = KMeans(n_clusters=2, n_init=10, random_state=seed + it).fit(Zi).labels_
        for a in range(n_s):
            for b in range(a + 1, n_s):
                key = (si[a], si[b])
                tot, same = cnt.get(key, (0, 0))
                cnt[key] = (tot + 1, same + int(lab[a] == lab[b]))
    return float(np.mean([1.0 - abs((s / t) - 0.5) * 2 for t, s in cnt.values()]))


def metrics_old(X_sub, seed):
    Z = np.asarray(X_sub, dtype=float)
    if Z.shape[1] < 3:
        return None
    Z = np.nan_to_num(Z, nan=0.0)
    Zn = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    S = PCA(n_components=3, random_state=SEED).fit_transform(Z)
    lab = KMeans(n_clusters=2, n_init=50, random_state=SEED).fit_predict(S[:, :3])
    sil = silhouette_score(S[:, :3], lab)
    b1 = GaussianMixture(1, covariance_type="diag", random_state=SEED).fit(S[:, :3]).bic(S[:, :3])
    b2 = GaussianMixture(2, covariance_type="diag", random_state=SEED).fit(S[:, :3]).bic(S[:, :3])
    return {"sil": float(sil), "bic_k1": float(b1), "bic_k2": float(b2),
            "bic_prefer": "k1" if b1 < b2 else "k2",
            "pac": pac_fixed(Zn, seed)}


def main():
    rows = []
    for cancer in ["GBM", "LGG"]:
        X = load_spladder_hivar(cancer)
        full_event = X.shape[1]
        tgt = min(TARGET[cancer], full_event)
        log(f"{cancer}: SplAdder 高变异事件池 {full_event}，下采样目标 {tgt}")

        full = metrics(X, SEED)
        log(f"  [全量 {full_event}] sil={full['sil']:.3f} PAC={full['pac']:.3f} "
            f"BIC偏好={full['bic_prefer']}")

        sims = []
        for rep in range(N_REP):
            cols = np.random.default_rng(SEED + rep).choice(
                X.columns, size=tgt, replace=False)
            m = metrics(X[list(cols)], SEED + rep)
            if m:
                sims.append(m)
                if (rep + 1) % 5 == 0:
                    log(f"    重复 {rep+1}/{N_REP}: sil={m['sil']:.3f} "
                        f"PAC={m['pac']:.3f} BIC={m['bic_prefer']}")

        arr_sil = np.array([s["sil"] for s in sims])
        arr_pac = np.array([s["pac"] for s in sims])
        k2 = sum(1 for s in sims if s["bic_prefer"] == "k2")
        ss = SPLICESEQ_VAL[cancer]
        rows.append({
            "cancer": cancer,
            "SplAdder_events_full": full_event,
            "SplAdder_sil_full": round(full["sil"], 3),
            "SplAdder_PAC_full": round(full["pac"], 3),
            "SplAdder_BIC_full": full["bic_prefer"],
            "downsampled_to": tgt,
            "down_sil_mean": round(float(arr_sil.mean()), 3),
            "down_sil_sd": round(float(arr_sil.std()), 3),
            "down_PAC_mean": round(float(arr_pac.mean()), 3),
            "down_BIC_k2_frac": round(k2 / len(sims), 2),
            "SpliceSeq_sil": ss["sil"],
            "SpliceSeq_PAC": ss["pac"],
            "SpliceSeq_BIC": ss["bic"],
        })
        log(f"  下采样 sil={arr_sil.mean():.3f}±{arr_sil.std():.3f} | "
            f"SpliceSeq sil={ss['sil']} | 差异解释度="
            f"{abs(arr_sil.mean()-full['sil'])/max(abs(ss['sil']-full['sil']),1e-9)*100:.0f}%")

    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/downsample_control.csv", index=False)
    log("\n" + R.to_string(index=False))
    log(f"\n输出 {OUT}/downsample_control.csv")


if __name__ == "__main__":
    main()
