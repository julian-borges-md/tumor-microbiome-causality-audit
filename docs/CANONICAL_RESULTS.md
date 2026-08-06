# RO-2026-008 | WP0 | Canonical results, 10-seed audit (wp0_core-1.0.0)

Seeds: 0-9. n = 360 samples, 150 taxa, 6 cancers, 3 batches, confounding 0.75.

| Metric | Cohort A (real signal) | Cohort B (ZERO signal) |
|---|---|---|
| Random-CV accuracy | +0.620 +- 0.025  [+0.578, +0.653] | +0.289 +- 0.019  [+0.261, +0.314] |
| T1 pass rate | 10/10 | 10/10 |
| T2 pass rate | 10/10 | 10/10 |
| T3 margin | +0.247 +- 0.024  [+0.194, +0.283] | -0.084 +- 0.020  [-0.111, -0.047] |
| T3 pass rate | 10/10 | 0/10 |
| T5a delta | +0.232 +- 0.042  [+0.163, +0.303] | -0.023 +- 0.024  [-0.067, +0.029] |
| T5a pass rate | 10/10 | 0/10 |
| **Verdict** | **10/10 PASS** | **0/10 PASS** |

Separation is complete: T3 margin and T5a delta ranges do not overlap
between cohorts across any of the 10 seeds.

HEADLINE: T1 and T2, the two significance tests most commonly reported in
this literature, passed on ZERO-SIGNAL data in 10 of 10 seeds.


## Detection floor, 10 seeds

| Effect (log units) | Sensitivity | Empirical false discovery |
|---|---|---|
| 0.10 | 0.01 +- 0.02 | 0.008 |
| 0.20 | 0.03 +- 0.03 | 0.007 |
| 0.40 | 0.24 +- 0.10 | 0.007 |
| 0.60 | 0.70 +- 0.14 | 0.013 |
| 0.80 | 0.98 +- 0.02 | 0.013 |
| 1.20 | 1.00 +- 0.00 | 0.011 |

Reliable floor (sensitivity >= 0.80) at approximately 0.8 log units, n = 360.
False discovery controlled throughout (0.007 to 0.013 against 0.05 nominal).


## Confounding sweep, 3 seeds

| Confounding | Cohort | Accuracy | T3 margin | T5a delta |
|---|---|---|---|---|
| 0.333 | A | 0.569 +- 0.042 | +0.406 | +0.211 +- 0.038 |
| 0.333 | B | 0.169 +- 0.022 | +0.006 | -0.041 +- 0.008 |
| 0.500 | A | 0.562 +- 0.044 | +0.308 | +0.169 +- 0.036 |
| 0.500 | B | 0.199 +- 0.013 | -0.055 | -0.048 +- 0.021 |
| 0.650 | A | 0.594 +- 0.050 | +0.265 | +0.174 +- 0.011 |
| 0.650 | B | 0.263 +- 0.033 | -0.067 | -0.014 +- 0.014 |
| 0.800 | A | 0.628 +- 0.035 | +0.225 | +0.261 +- 0.026 |
| 0.800 | B | 0.323 +- 0.025 | -0.080 | -0.029 +- 0.022 |
| 0.950 | A | 0.690 +- 0.016 | +0.219 | +0.319 +- 0.025 |
| 0.950 | B | 0.409 +- 0.005 | -0.061 | -0.036 +- 0.007 |

On ZERO-signal data, accuracy rises from 0.169 at no confounding to 0.409
at 0.95 confounding, against a no-information rate of 0.167. T5a delta
remains negative at every level.