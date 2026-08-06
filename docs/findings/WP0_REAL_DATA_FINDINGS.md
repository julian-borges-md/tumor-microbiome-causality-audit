# RO-2026-008 | WP0 | SIA_001
## First application to real data: The Cancer Microbiome Atlas

Dataset: TCMA (Dohlman et al., Cell Host & Microbe 2021), DOI 10.7924/r4bk1j35s.
WGS, tissue, sample level, decontaminated. Retrieved and analysed directly.
611 to 625 samples, 14,492 taxa, 5 TCGA projects, 3 sequencing centers.

---

## Finding 1 — TCMA survives the audit that synthetic zero-signal data failed

Cancer-type discrimination across five projects:

| Test | Value | Verdict |
|---|---|---|
| Random-CV accuracy | 0.840 (NIR 0.290, delta +0.550) | PASS, p = 8.6e-176 |
| T3 batch-only baseline | 0.617, margin +0.223 | PASS |
| T5a within-batch CV | 0.904 vs wbNIR 0.667, delta +0.237 | PASS |

Per-center within-batch deltas: +0.022, +0.250, +0.439. Heterogeneous; one
center carries almost no within-batch signal.

Interpretation: the decontaminated signal is real, not batch artifact.
Note the batch-only baseline reached 0.617 unaided, confirming substantial
center-outcome confounding is present in TCGA even after decontamination.

---

## Finding 2 — The signal is tissue of origin, not tumor state

Tumor vs patient-adjacent normal, WITHIN each cancer type (AUC, 0.5 = none):

| Cancer | Tumor n | Normal n | AUC | DA taxa (BH-FDR<0.05) |
|---|---|---|---|---|
| COAD | 125 | 21 | 0.506 | 0 / 916 |
| ESCA | 62 | 22 | 0.576 | 0 / 446 |
| STAD | 128 | 39 | 0.663 | 2 / 452 |
| HNSC | 157 | 22 | 0.693 | 51 / 913 |
| READ | 45 | 4 | not evaluable | - |

Cross-organ discrimination is strong; within-organ tumor-vs-normal
discrimination is weak to absent. Most of what reads as "tumor microbiome"
in pan-cancer analyses is anatomy.

---

## Finding 3 — Where taxa do differ, they are DEPLETED in tumor

| Cancer | Top differential taxa | Direction |
|---|---|---|
| HNSC | Rothia spp. (multiple strains), Actinomyces sp. ph3 | all depleted in tumor |
| STAD | Helicobacter suis, Helicobacter cetorum | depleted in tumor |

No taxon was significantly ENRICHED in tumor in any cancer type.
The HNSC count of 51 is inflated by strain-level redundancy (many Rothia
HMSC strains represent the same organism).

---

## Finding 4 (headline) — A proven causal carcinogen is DEPLETED at the tumor site

Helicobacter pylori in gastric cancer (STAD), an IARC Group 1 carcinogen and
the established cause of the disease:

| Metric | Tumor | Adjacent normal |
|---|---|---|
| Mean log-abundance | 0.158 | 1.200 |
| Prevalence | 0.48 | 0.54 |
| Mann-Whitney p | 0.022 | |

H. pylori is roughly 8-fold DEPLETED in the tumour relative to adjacent
normal tissue. This is biologically coherent: gastric carcinogenesis destroys
the mucosal niche (atrophy, intestinal metaplasia, loss of acid-secreting
cells) that H. pylori requires. The organism causes the cancer and is then
displaced by it.

### Consequence for causal inference
A known causal driver shows NEGATIVE cross-sectional association with the
tumour it caused. Therefore:
1. Cross-sectional tumour-vs-normal abundance CANNOT establish causation, and
   cannot even establish direction.
2. Absence of an organism from a tumour does not exclude it as the cause.
3. Enrichment of an organism in a tumour does not implicate it.
4. Reverse causation operates in BOTH directions: tumours can eliminate their
   own cause, not merely accumulate passengers.

This is empirical justification, from real data, for WP1's design decision to
anchor causal claims in germline instruments and temporal ordering rather than
intratumoural abundance.

---

## Positive control status

| Control | Result |
|---|---|
| Fusobacterium nucleatum in COAD | Directionally correct (tumour mean 0.191 vs normal 0.086, prevalence 0.58 vs 0.43) but NOT significant, p = 0.231 |
| H. pylori in STAD | Detected, significant, direction inverted (see Finding 4) |
| Escherichia coli | ZERO reads in every sample |
| Streptococcus gallolyticus | ZERO reads in every sample |

The Fusobacterium result is the expected direction at roughly 2.2-fold but
underpowered at n = 21 normals, consistent with the detection floor
established in Analysis 3. The complete absence of E. coli is a decontamination
artefact flag: E. coli is a classic reagent contaminant, and TCMA's filter
appears to have removed it globally, which would also remove any genuine
signal.

---

## Limitations (binding)

1. Normal-tissue sample sizes are 21 to 39. Underpowered by our own floor.
2. Analysis is unpaired; TCMA supports patient-matched pairing, not yet used.
3. Adjacent normal tissue is not healthy tissue (field cancerisation).
4. Single dataset, WGS only, five gastrointestinal and oropharyngeal cancers.
5. TCMA decontamination may over-filter (E. coli case).
6. Nulls here mean "no signal above the stated detection floor", not "no signal".
