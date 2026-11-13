#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""52 — 生成手稿表格（评审 P0-3：正文引用 Table 1/2 却无表）。

输出：
  manuscript/TABLES.md               Table 1（队列基线）、Table 2（CGGA 生存分层）
  manuscript/SUPPLEMENTARY_TABLES.md  Table S1（M0–M7×17 面板）、S2（PAC 模拟对照）、
                                      S3（残差化稳健性）
所有数字均来自本地实测 tsv，不做任何填充或估算。
"""
import os, json
import numpy as np
import pandas as pd

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma"
RI = os.path.join(ROOT, "ri_atlas")
OUTK = os.path.join(ROOT, "ri_atlas", "output_km")
OUTP = os.path.join(ROOT, "ri_atlas", "output_pancan")
MAN = os.path.join(ROOT, "manuscript")

MODELS = ["M0_raw", "M1_TotalExpr", "M2_plus_Neuronal", "M3_plus_Spliceosome",
          "M4_plus_Purity", "M5_plus_IDH", "M6_alt_Detection", "M7_hi_lo_layers"]
MLAB = {"M0_raw": "M0", "M1_TotalExpr": "M1", "M2_plus_Neuronal": "M2",
        "M3_plus_Spliceosome": "M3", "M4_plus_Purity": "M4", "M5_plus_IDH": "M5",
        "M6_alt_Detection": "M6", "M7_hi_lo_layers": "M7"}


def fmt_p(p):
    if not np.isfinite(p):
        return "—"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"


# ---------------------------------------------------------------- Table 1
def table1():
    rows = []
    # --- CGGA（本地临床，权威）---
    cg = pd.read_csv(os.path.join(RI, "output", "cgga_riaxis_scores.tsv"), sep="\t")
    for cohort, g in cg.groupby("cohort"):
        n = len(g)
        age = g["Age"].dropna()
        grade = g["Grade"].value_counts()
        rows.append({
            "cohort": cohort.replace("CGGA_", "CGGA-"),
            "n": n,
            "age": f"{age.median():.0f} ({age.quantile(.25):.0f}–{age.quantile(.75):.0f})",
            "male": f"{int((g['Gender'] == 'Male').sum())} ({(g['Gender'] == 'Male').mean() * 100:.0f}%)",
            "grade": f"{int(grade.get('II', 0))}/{int(grade.get('III', 0))}/{int(grade.get('IV', 0))}",
            "idh": f"{int((g['IDH'] == 'Mutant').sum())} ({(g['IDH'] == 'Mutant').mean() * 100:.0f}%)",
            "codel": f"{int((g['codel_1p19q'] == 'Codel').sum())} ({(g['codel_1p19q'] == 'Codel').mean() * 100:.0f}%)",
            "deaths": f"{int(g['event'].sum())} ({g['event'].mean() * 100:.0f}%)",
            "note": "expression-surrogate transfer of the IR axis (no PSI)",
        })

    # --- TCGA（严格限定为 IR 量化分析集；基线来自 cBioPortal pan-can-atlas 2018）---
    tb = "/tmp/tcga_baseline.json"
    # 分析集 barcode：fig3_scores.tsv 的样本（前 12 位 = patient barcode）
    sco = pd.read_csv(os.path.join(OUTK, "fig3_scores.tsv"), sep="\t")
    sco["key"] = sco["sample"].astype(str).str[:12]
    mol = pd.read_csv(os.path.join(OUTP, "glio_molecular.tsv"), sep="\t")
    for c in ("GBM", "LGG"):
        sub = sco[sco["cancer"] == c]
        keys = set(sub["key"])
        n_an = len(sub)                       # 样本数（手稿口径）
        n_pt = len(keys)                      # 病人数
        g = mol[mol["projectID"] == f"TCGA-{c}"]
        grade = (f"— / — / {n_an}" if c == "GBM" else "—")
        age = male = deaths = "—"
        if os.path.exists(tb):
            d = {k: v for k, v in json.load(open(tb)).get(c, {}).items() if k in keys}
            ages = [float(v.get("AGE")) for v in d.values()
                    if v.get("AGE") not in (None, "", "NA")]
            if ages:
                age = (f"{np.median(ages):.0f} "
                       f"({np.percentile(ages, 25):.0f}–{np.percentile(ages, 75):.0f})"
                       f" [n = {len(ages)}]")
            sex = [v.get("SEX") for v in d.values() if v.get("SEX")]
            if sex:
                nm = sum(s == "Male" for s in sex)
                male = f"{nm} of {len(sex)} ({nm / len(sex) * 100:.0f}%)"
            os_ = [v.get("OS_STATUS", "") for v in d.values()]
            os_ = [s for s in os_ if s]
            if os_:
                nd = sum("DECEASED" in s.upper() for s in os_)
                deaths = f"{nd} of {len(os_)} ({nd / len(os_) * 100:.0f}%)"
        idhm = int((g["IDH_status"] == "IDHmut").sum()) if len(g) else 0
        cd = int((g["codel_proxy"] == "codel_like").sum()) if len(g) else 0
        rows.append({
            "cohort": f"TCGA-{c}",
            "n": f"{n_an}" + (f" ({n_pt} patients)" if n_pt != n_an else ""),
            "age": age, "male": male,
            "grade": grade,
            "idh": f"{idhm} of {len(g)} annotated ({idhm / max(len(g), 1) * 100:.0f}%)" if len(g) else "—",
            "codel": f"{cd} of {len(g)} annotated ({cd / max(len(g), 1) * 100:.0f}%)" if len(g) else "—",
            "deaths": deaths,
            "note": "IR PSI analysed (primary cohort)",
        })

    df = pd.DataFrame(rows)
    order = ["TCGA-GBM", "TCGA-LGG", "CGGA-325", "CGGA-693"]
    df["_o"] = df["cohort"].map({k: i for i, k in enumerate(order)})
    df = df.sort_values("_o").drop(columns="_o")

    hdr = ("| Characteristic | " + " | ".join(df["cohort"]) + " |\n"
           "|" + "---|" * (len(df) + 1))
    def line(name, col):
        return f"| {name} | " + " | ".join(str(x) for x in df[col]) + " |"
    body = "\n".join([
        line("Cohort role", "note"),
        line("Tumours analysed, n", "n"),
        line("Median age at diagnosis, years (IQR)", "age"),
        line("Male sex, n (%)", "male"),
        line("WHO grade II / III / IV, n", "grade"),
        line("IDH-mutant, n (%)", "idh"),
        line("1p/19q codeletion (or proxy), n (%)", "codel"),
        line("Deaths during follow-up, n (%)", "deaths"),
    ])
    return ("**Table 1. Cohort characteristics.** "
            "TCGA cohorts are the primary IR-PSI analysis set; CGGA cohorts were used "
            "for the survival transfer only, using an expression-level surrogate of the "
            "IR axis rather than PSI. TCGA baseline age and sex were retrieved from "
            "cBioPortal (TCGA Pan-Cancer Atlas 2018); IDH status is from IDH1/IDH2 "
            "mutation and 1p/19q codeletion is proxied by co-occurring CIC/FUBP1 "
            "mutation, as described in the text. WHO grade is not reported per-tumour "
            "for TCGA-LGG in the retrieved attributes.\n\n" + hdr + "\n" + body)


# ---------------------------------------------------------------- Table 2
def table2():
    uni = pd.read_csv(os.path.join(RI, "output", "cgga_cox_univariate.tsv"), sep="\t")
    mul = pd.read_csv(os.path.join(RI, "output", "cgga_cox_multivariate.tsv"), sep="\t")
    lr = pd.read_csv(os.path.join(RI, "output", "cgga_logrank.tsv"), sep="\t")

    sub_cn = {"全队列": "Whole cohort", "WHO IV (GBM)": "WHO IV only",
              "WHO III": "WHO III only", "WHO II": "WHO II only",
              "IDH-wildtype": "IDH-wildtype only", "IDH-mutant": "IDH-mutant only",
              "IDH-wt WHO IV": "IDH-wildtype, WHO IV",
              "1p19q codel": "1p/19q codeleted", "1p19q non-codel": "1p/19q non-codeleted"}

    r = ["**Table 2. Survival transfer of the IR axis in the CGGA cohort "
         "(n = 970; expression-level surrogate, per SD increase).**",
         "Hazard ratios are univariable unless stated. The apparent whole-cohort "
         "association is not present in any stratum defined by grade, IDH or both, "
         "and does not survive multivariable adjustment.",
         "",
         "| Stratum | n | Events | HR per SD (95% CI) | p |",
         "|---|---:|---:|---|---:|"]
    for _, x in uni.iterrows():
        r.append(f"| {sub_cn.get(x['subset'], x['subset'])} | {int(x['n'])} | "
                 f"{int(x['events'])} | {x['HR_per_SD']:.3f} "
                 f"({x['lower']:.3f}–{x['upper']:.3f}) | {fmt_p(x['p'])} |")

    mv = mul[mul["term"] == "z"]
    mv_lab = {"全队列 校正分级/IDH/年龄/性别": "Whole cohort, adjusted for grade + IDH + age + sex",
              "全队列 +1p19q": "Whole cohort, additionally adjusted for 1p/19q",
              "IDH-wt WHO IV": "IDH-wildtype WHO IV, adjusted for age + sex"}
    r += ["", "**Multivariable adjustment (IR-axis score *z*):**", "",
          "| Model | n | Events | HR per SD (95% CI) | p |", "|---|---:|---:|---|---:|"]
    for _, x in mv.iterrows():
        r.append(f"| {mv_lab.get(x['subset'], x['subset'])} | {int(x['n'])} | "
                 f"{int(x['events'])} | {x['HR']:.3f} "
                 f"({x['lower']:.3f}–{x['upper']:.3f}) | {fmt_p(x['p'])} |")

    r += ["", "**Log-rank test of the tertile split:**", "",
          "| Stratum | n | χ² | p |", "|---|---:|---:|---:|"]
    for _, x in lr.iterrows():
        r.append(f"| {sub_cn.get(x['subset'], x['subset'])} | {int(x['n'])} | "
                 f"{x['chisq']:.2f} | {fmt_p(x['p'])} |")
    return "\n".join(r)


# ---------------------------------------------------------------- Table S1
def table_s1():
    bh = pd.read_csv(os.path.join(OUTP, "tme_partial_bh.tsv"), sep="\t")
    out = ["**Table S1. Rank partial correlations between IR burden and 17 "
           "microenvironment panels across the eight adjustment models.**",
           "M0 unadjusted; M1 +global expression; M2 +neuronal content; "
           "M3 +spliceosome core and proliferation; M4 +tumour content; M5 +IDH; "
           "M6 substitution of detection rate for global expression; "
           "M7 high/low expression layers + neuronal + spliceosome + purity + IDH. "
           "Entries are partial ρ; * q < 0.05 (BH within cohort and model).",
           ""]
    for cohort in ("GBM", "LGG"):
        g = bh[bh["cohort"] == cohort]
        panels = list(dict.fromkeys(g["panel"]))
        out.append(f"**{cohort}** (n = {int(g['n'].max())})")
        out.append("")
        out.append("| Panel | " + " | ".join(MLAB[m] for m in MODELS) + " |")
        out.append("|" + "---|" * (len(MODELS) + 1))
        piv = g.pivot_table(index="panel", columns="model", values="rho", aggfunc="first")
        qv = g.pivot_table(index="panel", columns="model", values="q", aggfunc="first")
        for p in panels:
            cells = []
            for m in MODELS:
                v = piv.loc[p, m] if (p in piv.index and m in piv.columns) else np.nan
                q = qv.loc[p, m] if (p in qv.index and m in qv.columns) else np.nan
                if not np.isfinite(v):
                    cells.append("—")
                else:
                    cells.append(f"{v:+.3f}" + ("*" if np.isfinite(q) and q < 0.05 else ""))
            out.append(f"| {p} | " + " | ".join(cells) + " |")
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------- Table S2
def table_s2():
    v = pd.read_csv(os.path.join(OUTK, "pac_algorithm_validation.tsv"), sep="\t")
    gt = {"A_two_clusters_strong": "Two clusters, well separated (n = 75 + 75)",
          "B_two_clusters_weak": "Two clusters, weakly separated (n = 75 + 75)",
          "C_single_gaussian_cloud": "Single Gaussian cloud (n = 150)",
          "D_continuous_gradient": "Continuous gradient, no gap (n = 150)"}
    r = ["**Table S2. Benchmark of the two PAC implementations against simulations "
         "with known ground truth.**",
         "Each scenario is a 150 × 20 matrix analysed identically to the tumour data "
         "(80% sample / 80% feature resampling, 100 iterations, k = 2). The legacy "
         "implementation divides the consensus matrix by the number of resampling "
         "iterations; the corrected implementation divides by the number of iterations "
         "in which each sample pair co-occurred. The legacy PAC cannot separate "
         "'two clusters' (0.497) from 'continuous gradient' (0.524); the corrected "
         "PAC does (0.000 vs 0.155). Silhouette separates the same scenarios and is "
         "unaffected by the error.",
         "",
         "| Scenario (ground truth) | PAC, legacy | PAC, corrected | Silhouette k = 2 |",
         "|---|---:|---:|---:|"]
    for _, x in v.iterrows():
        r.append(f"| {gt.get(x['scenario'], x['scenario'])} | {x['PAC_legacy']:.3f} | "
                 f"{x['PAC_monti']:.3f} | {x['silhouette_k2']:.3f} |")
    fd = json.load(open(os.path.join(OUTK, "cluster_diagnostics_fixed.json")))
    r += ["", "Observed values in the tumour data (corrected implementation):", "",
          "| Cohort | n | Silhouette k = 2 | PAC k = 2 | BIC prefers k = 1 |",
          "|---|---:|---:|---:|:--:|"]
    for c, d in fd.items():
        r.append(f"| {c} | {d['n_samples']} | {d['silhouette_k2']:.3f} | "
                 f"{d['PAC']['2']:.3f} | {'yes' if d['BIC_prefers_k1'] else 'no'} |")
    return "\n".join(r)


# ---------------------------------------------------------------- Table S3
def table_s3():
    d = pd.read_csv(os.path.join(OUTK, "residual_robustness.tsv"), sep="\t")
    conds = ["R0_raw", "R1_resid_TotalExpr", "R2_resid_full"]
    clab = {"R0_raw": "Unadjusted", "R1_resid_TotalExpr": "Residualised on global expression",
            "R2_resid_full": "Residualised on global expression + purity + neuronal"}
    metrics = [("PC1_var", "PC1 variance explained", 3),
               ("silhouette_k2", "Silhouette, k = 2", 3),
               ("PAC_k2", "PAC, k = 2", 3),
               ("rho_PC1_burden", "ρ(PC1, IR burden)", 3),
               ("rho_PC1_TotalExpr", "ρ(PC1, global expression)", 3),
               ("rho_PC1_Sm_core", "ρ(PC1, Sm core)", 3)]
    out = ["**Table S3. Robustness of the continuous axis to removal of global "
           "transcript output.**",
           "Every high-variability IR event was residualised on the listed covariates "
           "(rank-normalised) and the full diagnostic battery was repeated on the "
           "residual matrix. R1 removes global expression; R2 additionally removes "
           "tumour content and the neuronal/normal-brain signature. BIC preferred a "
           "single Gaussian component in every condition and both cohorts.",
           ""]
    for cancer in ("GBM", "LGG"):
        g = d[d["cancer"] == cancer]
        piv = g.pivot_table(index="metric", columns="condition", values="value",
                            aggfunc="first").apply(pd.to_numeric, errors="coerce")
        out.append(f"**{cancer}** (n = {int(g['n_samples'].max())}; "
                   f"{int(g['n_events'].max())} high-variability events)")
        out.append("")
        out.append("| Diagnostic | " + " | ".join(clab[c] for c in conds) + " |")
        out.append("|" + "---|" * (len(conds) + 1))
        for m, lab, nd in metrics:
            cells = []
            for c in conds:
                v = piv.loc[m, c] if (m in piv.index and c in piv.columns) else np.nan
                cells.append("—" if not np.isfinite(v) else f"{v:.{nd}f}")
            out.append(f"| {lab} | " + " | ".join(cells) + " |")
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    os.makedirs(MAN, exist_ok=True)
    main = "# Tables\n\n" + table1() + "\n\n" + table2() + "\n"
    open(os.path.join(MAN, "TABLES.md"), "w").write(main)
    sup = ("# Supplementary tables\n\n" + table_s1() + "\n\n"
           + table_s2() + "\n\n" + table_s3() + "\n")
    open(os.path.join(MAN, "SUPPLEMENTARY_TABLES.md"), "w").write(sup)
    print("DONE -> manuscript/TABLES.md, manuscript/SUPPLEMENTARY_TABLES.md")
