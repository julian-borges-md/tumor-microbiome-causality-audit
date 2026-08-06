# RO-2026-008 | WP0 | SIA_001
## Methods toolkit and calibration results (validated, pre-data)

Three controls, all executed and validated against ground truth before any
real data enters the pipeline. Seed 20260722. n = 480 to 660 samples,
6 cancer types, 3 sequencing batches, 150 to 250 taxa.

---

### Analysis 1 — Baseline audit validation
Two synthetic cohorts. A carries genuine taxon-cancer signal. B carries ZERO
biological signal; taxa respond only to batch.

| Test | Cohort A (real) | Cohort B (zero signal) | Discriminates |
|---|---|---|---|
| Random-CV accuracy | 0.677 | 0.332 | No (NIR = 0.167) |
| T1 no-information rate | PASS +0.511 | PASS +0.165 | No |
| T2 label permutation | PASS | PASS | No |
| T3 batch-only baseline | PASS +0.288 | FAIL -0.058 | Yes |
| T5a within-batch CV | PASS +0.265 | FAIL -0.022 | Yes |
| Verdict | PASS | FAIL | Correct |

T1 and T2, the two tests most commonly reported in this literature, both
passed on data containing no biology whatsoever.

Correction v1 to v2: leave-one-batch-out is invalid as a gate under strong
batch-outcome confounding (it removes an outcome class from training and
produced a false negative). Replaced by within-batch cross-validation.

---

### Analysis 2 — Confounding sensitivity sweep
How much apparent signal is manufactured from zero biology.

Cohort B (ZERO biological signal):

| Confounding | Random CV | Inflation over NIR | Batch-only | Within-batch delta |
|---|---|---|---|---|
| 0.333 (none) | 0.202 | +0.035 | 0.163 | -0.032 |
| 0.500 | 0.196 | +0.029 | 0.246 | -0.002 |
| 0.650 | 0.290 | +0.123 | 0.340 | -0.035 |
| 0.800 | 0.319 | +0.152 | 0.390 | +0.016 |
| 0.950 | 0.423 | +0.256 | 0.481 | +0.020 |

Cohort A (real signal):

| Confounding | Random CV | Batch-only | Within-batch delta |
|---|---|---|---|
| 0.333 | 0.640 | 0.163 | +0.297 |
| 0.500 | 0.656 | 0.246 | +0.228 |
| 0.650 | 0.627 | 0.340 | +0.192 |
| 0.800 | 0.608 | 0.390 | +0.284 |
| 0.950 | 0.752 | 0.481 | +0.368 |

Separation is complete across the full confounding range. Within-batch delta
never exceeds +0.020 on zero-signal data and never falls below +0.192 on real
signal.

---

### Analysis 3 — Detectability floor (n = 480)

| True effect (log units) | Sensitivity | Empirical false discovery | Within-batch delta | Verdict |
|---|---|---|---|---|
| 0.10 | 0.00 | 0.018 | -0.001 | below floor |
| 0.20 | 0.05 | 0.009 | +0.023 | below floor |
| 0.40 | 0.40 | 0.027 | +0.030 | partial |
| 0.60 | 0.90 | 0.018 | +0.081 | detectable |
| 0.80 | 1.00 | 0.027 | +0.143 | detectable |
| 1.20 | 1.00 | 0.027 | +0.270 | detectable |

Reliable detection floor at this sample size is an effect of approximately
0.6 log units. False discovery control holds throughout (0.9 to 2.7 percent
against a 5 percent nominal rate).

Known limitation: between effects of roughly 0.2 and 0.5, real signal exists
but falls below the audit's pass threshold. The audit is deliberately
conservative and will reject weak true signal. Any null must therefore be
reported as "no signal above the stated detection floor", never as
"no signal".

---

### Operating envelope derived
1. Any predictive claim requires T3 and T5a. T1 and T2 alone are worthless.
2. Signal survival must be reported against the detection floor.
3. Confounding strength must be measured and reported per cancer type.
