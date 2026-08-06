# RO-2026-008 | WP0 | Corrections to session-reported values

During development, simulation parameters were changed three times to fit
execution windows (taxa 500/300/250/150; permutations 200/40/25; trees
200/80/60/50). Values reported mid-session therefore came from different
configurations and are NOT mutually comparable. wp0_core-1.0.0 fixes one
canonical configuration and re-derives everything. The following corrections
apply.

## C1. Detection floor was overstated (single seed)

| Effect | Session value (1 seed, n=480) | Canonical (10 seeds, n=360) |
|---|---|---|
| 0.40 | sensitivity 0.40 | 0.24 +- 0.10 |
| 0.60 | sensitivity 0.90 | 0.70 +- 0.14 |
| 0.80 | sensitivity 1.00 | 0.98 +- 0.02 |

The reliable detection floor is an effect of approximately 0.8 log units at
n = 360, NOT 0.6 as stated during the session. Two causes: the canonical
configuration uses a smaller sample (360 vs 480), and the session value came
from a single favourable seed. Manuscript draft section 3.3 must be updated.

## C2. Synthetic point estimates superseded

All single-seed synthetic accuracies quoted during the session (for example
cohort A 0.677, cohort B 0.332) are superseded by the 10-seed canonical
values. Direction, separation and every PASS/FAIL verdict are unchanged.

## C3. Leave-one-batch-out invalidity (already documented)

LOBO was specified, found invalid under batch-outcome confounding, and
replaced by within-batch cross-validation. Retained as diagnostic only.

## C4. Redundancy collapsing rule (already documented)

First rule over-collapsed (41 organisms into one Actinobacteria clade);
correlation-gated rule at r>0.99 under-collapsed. Reported as a range,
not a point estimate.

## What did NOT change
All eleven real-data values verified exactly under wp0_tcma.py:
sample counts, cancer-type audit accuracy, batch-only baseline, within-batch
accuracy, all four paired-cohort sizes, and the H. pylori paired difference
and p-value.

---

# Appended 2026-08-05 | Audit of the retrieved package

## C5. Alistipes colorectal FinnGen effect size was misreported

The session narrative of 2026-07-29, and every record derived from it,
carried genus Alistipes in FinnGen colorectal cancer as OR 0.83, p 0.026.

| Source | OR | p (IVW) |
|---|---|---|
| Session narrative and checkpoint | 0.83 | 0.026 |
| MR_results_colorectal.tsv (canonical) | 0.95 displayed, 0.957 reconstructed | 0.578 |

No value resembling 0.83 appears in any column of that row, nor in any
related Rikenellaceae lineage row. The FinnGen colorectal test is null.

Consequence. The claim "nominally protective in all three independent tests"
is WITHDRAWN. The corrected claim is: protective direction in three of three,
nominal significance in TWO of three. The 100k colorectal and pancreatic
values are confirmed correct (0.933 p 0.048; 0.682 p 0.017).

The framing "lead, not a finding" is unaffected and remains correct.
See MR_CROSS_CONSISTENCY_FINDINGS.md section 5.

## C6. Directional agreement was tie-dependent as originally computed

The 61 percent cross-cohort agreement figure depended on the handling of
taxa whose odds ratio rounds to 1.00 in the FinnGen table.

| Tie rule | Agreement | Binomial p |
|---|---|---|
| Ties counted as disagreement | 56.9 percent | 0.054 |
| Ties dropped | 60.9 percent | 0.0027 |
| Sign recovered from log-CI midpoint (canonical) | 61.4 percent | 0.0011 |

The two naive rules straddle conventional significance. Sign recovery removes
the tie class entirely and confirms the original figure. The finding stands
and is now robust to the rule. Canonical value: 61.4 percent, 129 of 210.

## C7. Two reproducibility gaps closed, one residual

Closed: the cross-consistency analysis had no script and no findings record
(now wp1_cross_consistency.py and MR_CROSS_CONSISTENCY_FINDINGS.md); the
FinnGen colorectal table had no beta or standard error (now reconstructed
from the confidence interval, per-taxon QC in MR_cross_cohort_recon.tsv).

Residual: the FinnGen colorectal arm remains a RECONSTRUCTION, validated
against stored p-values but not re-derived from source. Definitive fix
requires re-running that MR with full-precision output, which needs an
OpenGWAS JWT. One taxon, phylum Cyanobacteria, fails reconstruction QC at
the p-value tail and is excluded; it was already diagnosed as a pleiotropic
artifact and its direction flipped on replication.
