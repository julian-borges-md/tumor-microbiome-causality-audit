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


---

# Appended 2026-08-07 | Figure regeneration audit

## C8. The 2.5x confounding ratio was attached to the wrong confounding level

README.md reported "2.5x chance at confounding strength 0.75". The sweep was
never run at 0.75. CANONICAL_RESULTS.md is correct: zero-signal accuracy rises
from 0.169 at 0.333 confounding to 0.409 at 0.95, against a no-information
rate of 0.167.

| Confounding | Zero-signal accuracy | Ratio to chance |
|---|---|---|
| 0.333 | 0.169 | 1.02x |
| 0.500 | 0.199 | 1.19x |
| 0.650 | 0.263 | 1.58x |
| 0.800 | 0.323 | 1.94x |
| 0.950 | 0.409 | **2.46x** |

Corrected claim: batch-outcome confounding alone drives accuracy to 2.5 times
chance on zero-signal data AT 0.95 CONFOUNDING, and to 1.9 times chance at
0.8. Figure 2 was always correct; only the README prose was wrong.

## C9. Within-tissue AUCs shift across scikit-learn versions

Figure 4 panel a was regenerated from the committed pipeline on scikit-learn
1.9.0. The cross-validated ROC AUCs move slightly from the values in the
original figure, produced on an earlier version.

| Cancer | Original figure | Regenerated | Samples |
|---|---|---|---|
| COAD | 0.52 | 0.502 | 125 v 21 |
| ESCA | 0.55 | 0.576 | 62 v 22 |
| STAD | 0.67 | 0.663 | 128 v 39 |
| HNSC | 0.68 | 0.693 | 157 v 22 |

Sample counts are identical. These AUCs were never entered into
CANONICAL_RESULTS.md and no claim depends on their second decimal. The claim
they support, weak within-tissue discrimination in all four evaluable cancers,
holds under both versions: every value lies between 0.50 and 0.70.

Figure 4 panel b re-derives exactly: paired H. pylori difference -0.9942
against the canonical -0.994, p = 4.46e-04, 39 matched pairs.

## C10. Committed source hardcoded the original container's path

Five scripts (wp0_tcma, wp0_tcma_real, wp0_paired, wp0_nesting,
wp0_tumor_vs_normal) contained `BASE = "/home/claude/tcma"`, an absolute path
that exists only in the container where they were written. They could not run
anywhere else. BASE now reads the TCMA_DIR environment variable and defaults
to ./tcma. Full real-data verification passes on a clean machine after a fresh
download: all checksums match and wp0_tcma.py reports VERIFICATION PASSED.


## C11. nest_sweep.py was not runnable as committed

The exploratory nesting sweep loaded its dependencies with
`exec(open('wp0_nesting.py').read().split('if __name__')[0])`, which required
being run from inside src/ and asserted nothing. Superseded by
wp0_nesting_sweep.py, which imports properly, emits
results/nesting_sweep.json, and asserts every value in
docs/findings/WP0_NESTING_FINDINGS.md. All values re-derive exactly:
STAD 8/3/3/2/2/2, HNSC 51/44/43/42/31/9, and the feature-space redundancy
table for all four cohorts. nest_sweep.py is retained unchanged as the
historical record.


## C12. MR instrument-selection parameters recovered from committed outputs

The MR arm had no committed code (limitation L2). Rather than assume the
parameters, they were recovered by testing candidate settings against the
committed n_snp counts for all 211 taxa.

The committed n_snp is the instrument count AFTER harmonisation with the
outcome GWAS, so any correct exposure-side setting must yield at least as many
instruments for EVERY taxon. A single deficit falsifies the setting.

| Parameters | Taxa with a deficit | Exact matches | Median selected/committed |
|---|---|---|---|
| p<1e-5, 500kb, palindromic dropped | **0** | 52/211 | **1.11** |
| p<1e-5, 1Mb, palindromic dropped | 0 | 52/211 | 1.10 |
| p<1e-5, 10Mb, palindromic dropped | 22 | 58/211 | 1.09 |
| p<1e-5, 500kb, palindromic kept | 0 | 17/211 | 1.23 |
| p<5e-6, any window | many | <10/211 | n/a |

The 10Mb window is falsified by its deficits. Keeping palindromic SNPs
inflates the ratio to 1.23, well above the ~10% loss expected from outcome
harmonisation. p<1e-5 with a 500kb to 1Mb distance window and palindromic
removal is the only setting consistent with the committed output, and it
matches the method as documented in MR_FINDINGS_colorectal.md.

Exposure data: MiBioGen MBG.allHits.p1e4.txt, 122,111 rows, 211 taxa,
sha256 37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e.
Publicly downloadable, no credentials.

### Status of L2

PARTIALLY closed. src/wp1_mr_pipeline.py commits the full method: instrument
selection, harmonisation, IVW, MR-Egger, weighted median, Cochran Q and
BH-FDR. Stage 1 (exposure) is validated against all 211 committed instrument
counts and runs today with no credentials. Stage 2 (outcome) is implemented
but NOT validated end to end, because the outcome summary statistics require
an OpenGWAS JWT (api.opengwas.io now returns 401 unauthenticated) or a direct
FinnGen / GWAS Catalog download that has not been performed. Stage 2 output
must not be used for any claim until that check is run. The module refuses to
run Stage 2 and says so.
