# Analysis code — *Global intron retention burden in glioma: a cross-cohort computational analysis of splicing structure, spliceosome coupling and clinical correlates*

This archive contains every script used to produce the results, tables and figures of the above
manuscript. It is deposited so that a reviewer can, in principle, re-derive every number in the
paper from the public input data.

**Licence:** CC-BY-4.0 · **Version:** v1.0 · **Contact:** mapenglaw@163.com (Peng Ma, Department of Neurosurgery, Taicang Loujiang New City Hospital)

---

## 1. Scope of what is archived

`code/` holds 60 scripts (44 Python, 16 R). They fall into eight analysis blocks, numbered so that
the intended execution order is the numeric order of the filename prefix:

| Prefix | Block | Key outputs |
|---:|---|---|
| `01` | IR burden definition and unsupervised structure (SplAdder/Kahles) | per-sample burden, PCA, k-means, GMM-BIC, silhouette, PAC |
| `02`–`07` | Clinical association, GSVA immune, survival, CGGA validation, snRNA localisation, figures | survival HRs, CGGA replication, single-cell maps |
| `08`–`14` | PSI–expression coupling, splicing-factor coupling, intron features, GSEA, mechanism figures | G1/G2, G3/G4, Fig 7, Fig 10 |
| `15`–`24` | Kahles-atlas lineage split: core metrics, confound control, survival, age positive control | GBM/LGG split analyses, Fig 5/6 |
| `30`–`43` | Pan-cancer scan (31 tumour types, 8,266 tumours), molecular subtypes, coupling | pan-cancer AUC, subtype comparison |
| `44`–`49` | Tumour-microenvironment partial coupling, scRNA paired dilution control | 17 immune panels, dilution permutation test |
| `50`–`54` | Residualisation robustness, PAC algorithm correction, table assembly, immune aggregate | Tables S1–S3 |
| `57` | **Replication with an independent quantification pipeline (TCGASpliceSeq)** | Table S4 |

**Superseded scripts.** `55_second_pipeline_spliceseq.py` and `56_downsample_control.py` are
exploratory and are retained for provenance only. They apply two conventions that were later found
to be inconsistent with the main pipeline and were corrected in `57`:

1. Burden was computed as the per-sample mean of the z-scored matrix. The manuscript definition —
   and the one used in `01` — is the **mean raw PSI across all IR events detected in ≥90% of
   tumours, with no SD filter**. Using the z-scored mean inflates ρ(PC1, burden) from 0.957 to 0.993.
2. Rank partial correlation used inverse-normal scores. The project convention (script `45`) is
   **pure ranks** (Pearson on ranks ≡ Spearman), with p from a t test on df = n − 2 − k. Inverse-normal
   scoring systematically attenuates coefficients.

**Two conventions that must not be changed** — altering either silently shifts reported values:

- **Burden** = mean raw PSI across all IR events with detection rate ≥ 0.90. No `SD ≥ 0.10` filter;
  that filter belongs to the clustering subset only. Verification: a recomputed burden must correlate
  with `output_pancan/tme_partial_samples.tsv::RI_burden` at Spearman ρ = 1.0000.
- **Rank partial correlation**: `rankdata` → OLS residualise both variables on the covariate design →
  Pearson on the residuals. With k = 0 covariates this equals Spearman exactly.

---

## 2. Input data (not included — all publicly available)

| Input | Source | Notes |
|---|---|---|
| Pan-cancer IR PSI matrices (SplAdder, hg19) | Kahles et al. 2018 pan-cancer splicing atlas | 123,390 IR events; GBM/LGG subsets cached as `ri_atlas/data/KM_IR_{GBM,LGG}.psi.gz` |
| Independent PSI quantification (hg19) | TCGASpliceSeq (Ryan et al., *Nucleic Acids Res* 2016;44:D1018-22) | used only by `57` |
| Gene-level expression, clinical and survival data | TCGA via UCSC Xena; CGGA (mRNAseq_325 / mRNAseq_693) | scripts `37`, `40` handle download and repair |
| Tumour purity, proliferation, neuronal/normal-brain content, IDH status | derived in scripts `17`, `34`, `45` | cached in `ri_atlas/output_pancan/tme_partial_samples.tsv` |
| Immune/microenvironment signature panels | published panels (17 in total), definitions in script `44` | — |
| Single-cell RNA-seq (paired dilution control) | public glioma scRNA-seq datasets, see script `47` | — |

> Raw data are not redistributed here because of their size and their own licence terms. Every input
> is obtainable from the cited public sources; `37_download_xena.py` and `40_finish_download.py`
> automate the TCGA portion.

---

## 3. Environment

Tested on macOS 15 (Apple Silicon) with:

- Python 3.13 — `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`
- R 4.6 — `ggplot2`, `dplyr`, `data.table`, `reshape2`, `scales`, `patchwork`, `ggpubr`,
  `survival`, `survminer`, `GSVA`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`,
  `ConsensusClusterPlus`

Minimal installation:

```bash
pip install numpy pandas scipy scikit-learn matplotlib
Rscript -e 'install.packages(c("ggplot2","dplyr","data.table","reshape2","scales",
  "patchwork","ggpubr","survival","survminer"), repos="https://cloud.r-project.org")'
```

**Random seeds.** All stochastic steps (k-means, Gaussian mixtures, PCA, subsampling) are pinned with
`random_state = SEED` where `SEED` is fixed per script; subsampling controls in `57` use
`np.random.default_rng(SEED)`. Reported values are therefore deterministic except for the 20-replicate
down-sampling control, where the manuscript reports mean ± SD across replicates.

---

## 4. Reproducing the headline numbers

| Manuscript claim | Script | Command |
|---|---|---|
| Burden definition, PC1, clustering diagnostics | `01_ri_burden_cluster.py` | `python 01_ri_burden_cluster.py` |
| Pan-cancer continuous-axis scan (8,266 tumours) | `30`–`36` | run in numeric order |
| Residualisation on global expression | `50_residual_robustness.py` | `python 50_residual_robustness.py` |
| PAC estimator correction (Table S2) | `51_pac_algorithm_fix.py` | `python 51_pac_algorithm_fix.py` |
| Microenvironment partial coupling (Table S1) | `45_tme_partial_coupling.py`, `54_immune_aggregate.py` | in order |
| **Independent-pipeline replication (Table S4)** | `57_cross_pipeline_replication.py` | `python 57_cross_pipeline_replication.py` |

Script `57` is the one most likely to be run first by a reviewer. It is self-contained, takes roughly
two minutes on a laptop, and writes `output/cross_pipeline_replication/{T1_results,
T2_cross_pipeline_results,T3_sm_results}.csv`. Its expected values are reproduced in Table S4 of the
manuscript; the burden recomputation check prints `Spearman ρ = 1.0000` against the archived main-pipeline
burden as a built-in assertion of definition consistency.

---

## 5. Known limitations of the archived code

- Paths inside the scripts are absolute and point to the original working machine. Before running,
  replace the `ROOT` / `DATA` constants at the head of each script with your own paths.
- Scripts `37`, `40`, `41`, `42` perform network downloads and repair partially transferred files;
  they are idempotent but slow and are not needed if data are already in place.
- `55` and `56` are superseded; see section 1.
