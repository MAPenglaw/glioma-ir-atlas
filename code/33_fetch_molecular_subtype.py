#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""33 — 获取 GBM/LGG 的分子分型（B-1）。

本地 TCGA_GBMLGG_clinical.tsv 缺 IDH / 1p19q 状态，这里从 cBioPortal 补。

数据源：cBioPortal REST（TCGA PanCanAtlas 2018）
  mutations profile，Entrez：
    IDH1 3417 / IDH2 3418 / TP53 7157 / ATRX 546 / EGFR 1956 /
    CIC 23152 / FUBP1 8880 / NF1 4763 / PTEN 5728

分型规则（遵循 WHO 2021 / Ceccarelli 2016 的公开做法）：
  IDHmut  = IDH1 或 IDH2 有非同义突变
  IDHwt   = 两者均无
  codel 代理 = IDHmut 且（CIC 或 FUBP1 突变）
      —— CIC(19q13.2) 与 FUBP1(1p31.1) 突变几乎只见于 IDHmut 伴 1p/19q 共缺失的
         少突胶质瘤，是文献常用的共缺失代理标记。cBioPortal 的 armlevel_cna
         接口在本环境 404，无法直接取臂水平缺失，故用此代理并在文中说明。

注意：只保留经 OncoKB/数据库过滤后的体细胞突变记录；cBioPortal 的
PanCanAtlas mutations 已做过常见胚系变异过滤。

输出：output_pancan/glio_molecular.tsv
  sample_id(16位) / patient(12位) / project / IDH_status / codel_proxy /
  各基因突变 0-1
"""
import os
import json
import urllib.request

ROOT = "/Users/peng/WorkBuddy/2026-09-21-10-00-22/snrpc_glioma/ri_atlas"
OUT = os.path.join(ROOT, "output_pancan")
CACHE = os.path.join(ROOT, "data", "pancan_clinical")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

GENES = {
    3417: "IDH1", 3418: "IDH2", 7157: "TP53", 546: "ATRX",
    1956: "EGFR", 23152: "CIC", 8880: "FUBP1", 4763: "NF1",
    5728: "PTEN",
}
STUDIES = {
    "LGG": "lgg_tcga_pan_can_atlas_2018",
    "GBM": "gbm_tcga_pan_can_atlas_2018",
}


# 本机有 HTTP 代理，会让 cBioPortal 请求失败（此前 Ensembl 同样中招），显式绕过
urllib.request.install_opener(
    urllib.request.build_opener(urllib.request.ProxyHandler({})))


def fetch(study, cache_name):
    path = os.path.join(CACHE, cache_name)
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"[cache] {cache_name} {os.path.getsize(path)} B", flush=True)
        return json.load(open(path))
    url = (f"https://www.cbioportal.org/api/molecular-profiles/"
           f"{study}_mutations/mutations/fetch?projection=SUMMARY")
    body = json.dumps({
        "sampleListId": f"{study}_sequenced",
        "entrezGeneIds": list(GENES.keys()),
    }).encode()
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/json"})
    print(f"[get] {study} mutations ...", flush=True)
    with urllib.request.urlopen(req, timeout=300) as r:
        data = json.loads(r.read().decode())
    json.dump(data, open(path, "w"))
    print(f"[get] {study} {len(data)} 条突变", flush=True)
    return data


def main():
    rows = {}
    for proj, study in STUDIES.items():
        data = fetch(study, f"{proj}_mut9.json")
        rec = {}
        for m in data:
            sid = m.get("sampleId") or ""
            pid = m.get("patientId") or (sid[:12] if sid else "")
            gid = m.get("entrezGeneId")
            sym = GENES.get(gid)
            if sym is None or not pid:
                continue
            # 过滤掉 Silent / Intron / 3'UTR 等非蛋白改变突变
            vtype = (m.get("variantType") or "").upper()
            if vtype in ("SILENT", "INTRON", "3'UTR", "5'UTR", "FLANK"):
                continue
            rec.setdefault(pid, {g: 0 for g in GENES.values()})
            rec[pid][sym] = 1
        print(f"[{proj}] 有突变记录的患者 {len(rec)}", flush=True)
        for pid, g in rec.items():
            rows[pid] = {"patient": pid, "project": proj, **g}

    import pandas as pd
    df = pd.DataFrame.from_records(list(rows.values()))
    for c in GENES.values():
        if c not in df.columns:
            df[c] = 0
    df = df.fillna(0)

    df["IDH_status"] = df.apply(
        lambda r: "IDHmut" if (r["IDH1"] == 1 or r["IDH2"] == 1) else "IDHwt",
        axis=1)
    df["codel_proxy"] = df.apply(
        lambda r: ("codel_like" if (r["IDH_status"] == "IDHmut" and
                                    (r["CIC"] == 1 or r["FUBP1"] == 1))
                   else ("noncodel_like" if r["IDH_status"] == "IDHmut"
                         else "NA")), axis=1)
    df["mol3"] = df.apply(
        lambda r: (r["codel_proxy"] if r["IDH_status"] == "IDHmut"
                   else "IDHwt"), axis=1)

    # 与 B 线样本对齐（16 位 barcode → 12 位患者）
    samp = pd.read_csv(os.path.join(OUT, "mat_hq_samples.tsv"), sep="\t")
    glio = samp[samp.projectID.isin(["TCGA-GBM", "TCGA-LGG"])].copy()
    glio["patient"] = glio["sample_id"].str[:12]
    mg = glio.merge(df, on="patient", how="left", suffixes=("", "_m"))

    print(f"\n=== 对齐结果 ===", flush=True)
    print(f"胶质瘤样本 {len(glio)}，其中匹配到突变记录 "
          f"{int(mg['IDH_status'].notna().sum())}", flush=True)
    print("\n各队列 IDH 状态：", flush=True)
    print(mg.groupby(["projectID", "IDH_status"]).size().to_string(), flush=True)
    print("\n三分类：", flush=True)
    print(mg.groupby(["projectID", "mol3"]).size().to_string(), flush=True)

    out = mg[["sample_id", "patient", "projectID", "IDH_status",
              "codel_proxy", "mol3"] + list(GENES.values())]
    out.to_csv(os.path.join(OUT, "glio_molecular.tsv"), sep="\t", index=False)
    df.to_csv(os.path.join(OUT, "glio_molecular_bypatient.tsv"),
              sep="\t", index=False)
    print(f"\n[sav] glio_molecular.tsv  {len(out)} 行", flush=True)


if __name__ == "__main__":
    main()
