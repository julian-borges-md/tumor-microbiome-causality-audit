[![verify](https://github.com/julian-borges-md/tumor-microbiome-causality-audit/actions/workflows/verify.yml/badge.svg)](https://github.com/julian-borges-md/tumor-microbiome-causality-audit/actions/workflows/verify.yml)
[![Data DOI](https://img.shields.io/badge/Dataset-TCMA_10.7924%2Fr4bk1j35s-blue)](https://doi.org/10.7924/r4bk1j35s)
[![Python](https://img.shields.io/badge/Code-Python-yellowgreen)](#)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-assertion--checked-brightgreen)](docs/MANIFEST.md)
[![Presented](https://img.shields.io/badge/Presented-BU_Health_Data_Science_%26_AI_Showcase_2026-cc0000)](https://sites.bu.edu/healthdatascience/)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

---

# Cross-Sectional Intratumoral Microbial Abundance Cannot Establish Causation

### A calibration study using an established carcinogen as a natural control

**Author:** Julian Borges, MD, MS
**Presented at:** BU Health Data Science & AI Showcase, 15 September 2026, Hiebert Lounge, Boston University Medical Campus

---

## Project Overview

Decontamination pipelines have made real progress on one question in tumor microbiome research: **is a detected taxon actually present?** They do not address a second question: **if it is present, is it causal, and in which direction does the association run?**

This project answers the second question by refusing to assume the instruments work. It builds ground-truth data where the answer is known in advance, tests whether the field's standard significance tests can detect the *absence* of signal, and only then touches real patients.

The headline result is uncomfortable. **The two significance tests most commonly reported in this literature both pass on data constructed to contain no biological signal, in ten of ten random seeds.** Applied to a decontaminated reference cohort, an IARC Group 1 carcinogen appears on the wrong side of the tumor-versus-normal comparison.

This is a direct methodological successor to [`breast-cancer-ai-misclassification-`](https://github.com/julian-borges-md/breast-cancer-ai-misclassification-) (*JAMIA Open*, 2026), which audited shortcut learning in genomic subtyping. Same failure mode, different substrate, stronger control.

---

## Findings

| # | Finding | Value |
|---|---|---|
| F1 | Standard significance tests pass on zero-signal data | **10 of 10 seeds.** Only a confounder baseline (T3) and within-batch cross-validation (T5a) discriminate |
| F2 | Batch-outcome confounding alone manufactures accuracy from nothing | **2.46x chance** at 0.95 confounding (accuracy 0.409 vs 0.167 no-information rate); 1.94x at 0.8 |
| F3 | *H. pylori*, the accepted cause of gastric adenocarcinoma, is **depleted** at the tumor site | **~8-fold**, 39 patient-matched pairs, difference &minus;0.99 log units, p = 4.5e&minus;4 |
| F4 | Reported discovery counts depend on the taxonomic redundancy rule | Up to **5.7-fold** variation on identical data at identical FDR |
| F5 | Detection floor | **~0.8 log units** at n = 360. A null below this is uninterpretable |
| F6 | MR screen of 211 gut taxa against colorectal cancer | **Zero taxa survive FDR** in a 100,204-case cohort. Every FinnGen lead failed replication |
| F7 | Cross-cohort directional agreement | **61.4%** (129/210, binomial p = 0.0011) against 50% under noise. Real but diffuse causal signal |

**What this does not show.** It does not show that published tumor microbiome findings are wrong. It shows that the tests as reported cannot distinguish signal from artifact. Some of that work may be entirely correct.

---

## Objectives

1. Build a ground-truth simulator emitting a signal cohort and a zero-signal cohort from the same code path.
2. Determine whether the audit battery can *fail* before running it on real data.
3. Quantify the detection floor rather than reporting an uncalibrated null.
4. Apply the validated battery to a decontaminated reference cohort.
5. Test whether cross-sectional abundance recovers the correct causal direction for an organism whose causal role is independently established.
6. Substitute germline anchoring where cross-sectional inference fails, with independent replication.

---

## Methods & Tools

| Layer | Tool |
|---|---|
| Language | Python 3.11+ |
| Models | scikit-learn random forests. Deliberately unexciting: if the finding depended on the model, the finding would be about the model |
| Numerics | numpy, scipy |
| Figures | matplotlib, Okabe-Ito colorblind-safe palette, 300 dpi |
| Tumor data | [The Cancer Microbiome Atlas](https://doi.org/10.7924/r4bk1j35s) (Dohlman 2021). 611 samples, 14,492 taxa, five TCGA projects, three sequencing centers |
| Exposure data | MiBioGen, 211 taxa, 18,340 individuals |
| Outcome data | FinnGen R12; Fernandez-Rozadilla 2023 (GCST90129505, 100,204 colorectal cancer cases) |
| Access | OpenGWAS REST API (JWT authentication) |
| Runtime | Single machine, no GPU, minutes not hours |

**Audit battery.** Five tests. T1 comparison against a no-information rate and T2 label permutation are what the field reports. T3 confounder baseline, T4 and T5a within-batch cross-validation were added here. Every test returns pass or fail, not a p-value to be interpreted after the fact.

**Causal arm.** Instrument selection, harmonization, inverse variance weighting with MR-Egger, weighted median, Cochran Q, Benjamini-Hochberg correction, and leave-one-out sensitivity.

---

## Reproduce

```bash
git clone https://github.com/julian-borges-md/tumor-microbiome-causality-audit.git
cd tumor-microbiome-causality-audit
pip install -r requirements.txt
make verify
```

`make verify` re-derives every asserted number in the manuscript from the committed result tables. It requires no network access and no raw data download. **The module exits non-zero and prints `DRIFT` to stderr if any value disagrees with the paper.** This is the check that runs in CI, so a red badge above means the repository no longer agrees with the manuscript.

Expected output:

```
Reconstruction QC: 210/211 pass, 1 fail
Directional agreement (sign-recovered): 129/210 = 61.4%, binomial p = 0.001131
Consistent-direction taxa across all three tests: 82/210 (expected under independence: 52.5)
All asserted values re-derived within tolerance.
```

Verified reproducing identically across Python 3.11/numpy 1.26/scipy 1.11 and Python 3.14/numpy 2.4/scipy 1.17.

The WP0 simulation and real-data targets require the TCMA download described in [`docs/RUNBOOK.md`](docs/RUNBOOK.md). Every claim maps to a command and an expected value in [`docs/MANIFEST.md`](docs/MANIFEST.md).

---

## Repository Contents

```
figures/make_figures.py       regenerates Figures 1, 2, 3, 6 from results/
figures/make_figures_real.py  regenerates Figures 4, 5 from TCMA (needs TCMA_DIR)
  wp0_nesting_sweep.py     canonical taxonomic redundancy sweep
src/                       analysis modules
  wp0_core.py              ground-truth simulator, signal and zero-signal cohorts
  wp0_confound_sweep.py    accuracy against batch-outcome confounding strength
  wp0_detection_floor.py   sensitivity by effect size, ten seeds
  wp0_tcma_real.py         audit battery applied to real decontaminated data
  wp0_paired.py            patient-matched tumor versus adjacent-normal analysis
  wp0_nesting.py           taxonomic redundancy sweep
  wp1_cross_consistency.py cross-cohort consistency, assertion-checked
results/                   committed result tables and JSON
figures/                   Figures 1-7, 300 dpi
figures/poster/            print poster, social cards, and their generators
docs/                      runbook, manifest, corrections log, findings records
```

---

## Reproducibility Standard

| Guarantee | Implementation |
|---|---|
| Determinism | Ten fixed seeds, `random_state` pinned |
| Provenance | SHA256 checksums on every input ([`docs/DATA_CHECKSUMS.txt`](docs/DATA_CHECKSUMS.txt)) |
| Claim traceability | Every claim maps to a command and an expected value ([`docs/MANIFEST.md`](docs/MANIFEST.md)) |
| Drift detection | Expected values compiled into modules as assertions. Failure exits non-zero |
| Error disclosure | Append-only corrections log ([`docs/CORRECTIONS.md`](docs/CORRECTIONS.md)) |

### Known limitations, disclosed

This repository documents its own defects rather than omitting them.

| ID | Limitation |
|---|---|
| L1 | **Closed.** All seven figures regenerate from source. Figures 1, 2, 3, 6 via `figures/make_figures.py`; Figures 4, 5 via `figures/make_figures_real.py`; Figure 7 via `src/wp1_cross_consistency.py`. Every module is assertion-checked and exits non-zero on drift |
| L2 | The Mendelian randomization pipeline is not committed. Result tables and plots only. Scheduled, not done |
| L3 | The FinnGen colorectal arm is **reconstructed** from rounded confidence intervals, validated against stored p-values, not re-derived from source. Definitive fix requires re-running with full-precision output |
| L4 | Permutation count is 25, not 1000, chosen for runtime |
| L5 | Classifier variance is not sampled (`random_state` fixed), so reported dispersion understates total uncertainty |
| L6 | Phylum Cyanobacteria is excluded from the agreement statistic on reconstruction QC failure at the p-value tail |
| L7 | MR limitations carry through: liberal p < 1e&minus;5 instrument threshold, distance-based clumping, palindromic SNPs dropped, partial ancestry mismatch, no colocalization, no reverse-direction MR |

**A correction is on record.** An effect size in the working notes (OR 0.83, p 0.026) did not exist in the saved data (OR 0.95, p 0.578). Root causes were an analysis that ran inline and was never persisted, and a result table that stored a formatted display string instead of the underlying estimate. Both are fixed. See [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) C5 through C7.

---

## Data Availability

| Source | Access | Identifier |
|---|---|---|
| The Cancer Microbiome Atlas | Open | [10.7924/r4bk1j35s](https://doi.org/10.7924/r4bk1j35s) |
| MiBioGen | Open | Kurilshikov 2021 |
| FinnGen R12 | Open, via OpenGWAS | finngen_R12_* |
| Fernandez-Rozadilla 2023 | Open, via GWAS Catalog | GCST90129505 |

No controlled-access or patient-identifiable data is used or redistributed here.

---

## Citation

> Borges, Julian. *Cross-sectional intratumoral microbial abundance cannot establish causation: a calibration study using an established carcinogen as a natural control.* 2026. https://github.com/julian-borges-md/tumor-microbiome-causality-audit

Manuscript under preparation. Zenodo archive and DOI on acceptance. Machine-readable metadata in [`CITATION.cff`](CITATION.cff).

---

## Author

**Julian Borges, MD, MS**
Department of Computer Science, Boston University
ORCID [0009-0001-9929-3135](https://orcid.org/0009-0001-9929-3135)
jyborges@bu.edu

---

## Keywords

`tumor microbiome` • `causal inference` • `shortcut learning` • `negative control` • `Mendelian randomization` • `reproducibility` • `detection floor` • `confounding` • `TCGA` • `health AI safety`

---

## License

Code is licensed under the [MIT License](LICENSE). Figures and documentation are released under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). Upstream datasets retain their original licenses.

---

<div align="center">

**Frontier Translational Research Lab**

Department of Computer Science · Boston University · Harvard Medical School GCSRT Alumni

[![Lab Website](https://img.shields.io/badge/Lab-frontier--lab-002244?style=flat-square)](https://julian-borges-md.github.io/frontier-lab/)
[![BU CS](https://img.shields.io/badge/BU-Computer_Science-cc0000?style=flat-square)](https://www.bu.edu/cs/)
[![HMS](https://img.shields.io/badge/HMS-GCSRT_Alumni-a51c30?style=flat-square)](https://ghsm.hms.harvard.edu/education/global-clinical-scholars-research-training)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0001--9929--3135-a6ce39?style=flat-square&logo=orcid&logoColor=white)](https://orcid.org/0009-0001-9929-3135)
[![CV](https://img.shields.io/badge/Academic_CV-research--profile-4f46e5?style=flat-square)](https://julian-borges-md.github.io/research-profile/)

*Julian Borges, MD, MS · jyborges@bu.edu*

</div>
