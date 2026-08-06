# RO-2026-008 | WP0 | Reproducibility Manifest
wp0_core-1.0.0

Every claim maps to a command and an expected value. Nothing reported is
un-derivable from this package.

## Package contents

| File | Role |
|---|---|
| RUNBOOK.md | Exact reproduction procedure, start to finish |
| ENVIRONMENT.txt | Interpreter and library versions |
| DATA_CHECKSUMS.txt | SHA256 of all input data, verify before running |
| wp0_core.py | Canonical synthetic analyses (audit, sweep, floor) |
| wp0_tcma.py | Canonical real-data analyses WITH verification assertions |
| CANONICAL_RESULTS.md | All multi-seed results, superseding session values |
| CORRECTIONS.md | Every value that changed and why |
| audit_seeds_0_4.json, audit_seeds_5_9.json | Raw 10-seed audit output |
| floor_10seeds.json, sweep_3seeds.json | Raw sweep and floor output |

## Claim to evidence map

| # | Claim | Command | Expected |
|---|---|---|---|
| 1 | Conventional tests pass on zero-signal data | `wp0_core.py audit --seeds 0..9` | T1 and T2 pass 10/10 in cohort B |
| 2 | Confounder baseline discriminates | same | T3 pass 10/10 in A, 0/10 in B |
| 3 | Within-batch CV discriminates | same | T5a pass 10/10 in A, 0/10 in B |
| 4 | Separation is complete | same | T3 and T5a ranges non-overlapping across all seeds |
| 5 | Confounding manufactures signal monotonically | `wp0_core.py sweep --seeds 0,1,2` | Cohort B accuracy 0.169 to 0.409 as confounding rises |
| 6 | Detection floor approx 0.8 log units | `wp0_core.py floor --seeds 0..9` | sensitivity 0.98 at 0.80, 0.70 at 0.60 |
| 7 | TCMA passes the audit | `wp0_tcma.py` | acc 0.840, batch-only 0.617, within-batch 0.904 |
| 8 | H. pylori depleted in gastric tumour | `wp0_tcma.py` | paired diff -0.994, p 4.5e-04, 39 pairs |
| 9 | Discovery counts are rule-dependent | `nest_sweep.py` (exploratory, not yet canonicalised) | STAD 8 to 2, HNSC 51 to 9 |

## Verification status

| Category | Status |
|---|---|
| Real-data values | ALL 11 VERIFIED exactly by assertion (wp0_tcma.py exits 0) |
| Synthetic values | RE-DERIVED under canonical config; session point estimates superseded |
| Detection floor | CORRECTED, was overstated from a single seed (see CORRECTIONS.md C1) |
| Nesting analysis | NOT yet canonicalised; still uses exploratory exec-based import |

## Known reproducibility gaps (honest)

1. **Nesting analysis (Analysis 7) is not in the canonical package.** It still
   relies on `exec(open(...).read().split(...))` string-splitting, which is
   fragile. It must be ported into wp0_tcma.py before submission.
2. **Classifier variance is not sampled.** random_state is fixed at 0
   throughout, so reported dispersion reflects simulation variance only and
   understates total uncertainty.
3. **Single dataset.** All real-data claims rest on TCMA. No independent
   cohort has been analysed.
4. **Permutation count is 25**, chosen for runtime. Adequate for a 95th
   percentile threshold but coarse; production should use 1000.
5. **Simulation parameters are illustrative**, not calibrated to any real
   assay's noise structure. They must not be read as a power calculation.
6. **Figures not generated.** No figure code exists yet.

---

# Appended 2026-08-05 | WP1 cross-consistency module

| Claim | Command | Expected |
|---|---|---|
| Directional agreement 61.4 percent, 129 of 210, binomial p 0.0011 | `python3 wp1_cross_consistency.py` | asserted, exit 0 |
| Reconstruction QC 210 of 211 pass, 1 fail | same | asserted |
| Alistipes CRC FinnGen OR 0.95 displayed, p 0.578 | same | asserted |
| Alistipes CRC 100k OR 0.9334, p 0.0478 | same | asserted |
| Alistipes pancreatic OR 0.6823, p 0.0168 | same | asserted |
| Triangulation 82 of 210 consistent direction | same | reported, descriptive only |

Module exits non-zero and prints DRIFT lines on any failure. Verified
deterministic across repeated runs.

## Reproducibility gap status update

Gap 6 in the list above, "Figures not generated. No figure code exists yet,"
is partially closed: figure code for Figure 7 is in wp1_cross_consistency.py.
Figures 1 to 6 remain without committed generation code.

New gap: the FinnGen colorectal MR arm is reconstructed from rounded
confidence intervals, not re-derived. See CORRECTIONS.md C7.
