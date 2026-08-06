# RO-2026-008 | WP1 | Cross-cohort consistency and triangulation

Module: `wp1_cross_consistency-1.0.0`
Derived from: `MR_results_colorectal.tsv`, `MR_results_colorectal_100k.tsv`, `MR_results_pancreatic.tsv`
Reproduce with: `python3 wp1_cross_consistency.py` (exit 0 required)

This document canonicalises an analysis that was executed inline on 2026-07-29
and never persisted. It supersedes every value reported for this analysis in
the session narrative and in the 2026-07-29 checkpoint.

---

## 1. Why this analysis needed rebuilding

Three defects were found on audit of the retrieved package.

| ID | Defect | Resolution |
|---|---|---|
| D2 | The reported 61 percent directional agreement depended on how taxa with an odds ratio rounding to 1.00 were handled. Counting them as disagreements gave 56.9 percent at p 0.054; dropping them gave 60.9 percent at p 0.0027. The two rules straddle conventional significance | Sign recovered for every taxon from the log confidence interval midpoint, removing the tie class entirely. Both naive rules retained and reported as sensitivity |
| D3 | No script and no findings record existed for the analysis. It failed the standard set in MANIFEST.md, that every claim maps to a command and an expected value | This module and this document |
| D4 | `MR_results_colorectal.tsv` stores the effect as a rounded display string with no beta and no standard error, so the FinnGen colorectal arm was not re-analysable | Beta and standard error reconstructed from the confidence interval and validated against the stored p-value per taxon. Reconstruction written to `MR_cross_cohort_recon.tsv` with a QC column |

---

## 2. Reconstruction method and its validation

For the FinnGen colorectal table the log-scale effect is recovered as the
midpoint of the log confidence interval rather than as the log of the
displayed odds ratio:

    b_hat  = (ln LCI + ln UCI) / 2
    se_hat = (ln UCI - ln LCI) / (2 x 1.959964)

The midpoint carries more information than the two decimal odds ratio because
it is built from two independently rounded bounds. A taxon displayed as
OR 1.00 with an asymmetric interval still yields a signed effect.

Each reconstruction is checked by recomputing the two sided p-value and
comparing it to the stored value on the log10 scale, tolerance 0.15.

| Result | Value |
|---|---|
| Taxa reconstructed | 211 |
| Passed QC | 210 |
| Failed QC | 1 |

The single failure is **phylum Cyanobacteria**, which carries the smallest
p-value in the table at 1.16e-4. Confidence interval rounding loses precision
in the tail, so its reconstructed p of 6.8e-5 falls outside tolerance. It is
excluded rather than used. This costs nothing analytically: that taxon was
already diagnosed as a pleiotropic artifact through the NOS2 instrument
rs2314810, and its direction flipped on replication.

---

## 3. Result 1: directional agreement between two independent colorectal cohorts

210 taxa tested in both FinnGen R12 (11,790 cases) and Fernandez-Rozadilla
2023 (100,204 cases).

| Tie rule | Agreement | Binomial p versus 0.5 |
|---|---|---|
| **Sign recovered (canonical)** | **129 of 210 = 61.4 percent** | **0.0011** |
| Naive, ties dropped | 61.2 percent | 0.0021 |
| Naive, ties counted as disagreement | 57.1 percent | 0.045 |

**The finding is confirmed and is now robust to the tie rule.** All three rules
exceed the 50 percent noise expectation and all three reach nominal
significance. The canonical value is 61.4 percent.

Interpretation is unchanged from the original reading and remains the correct
one. This is not evidence of a causal driver. It is evidence that the causal
signal from the gut microbiome on colorectal cancer is real but diffuse,
spread thinly across many taxa, with no single organism strong enough to
survive multiple testing correction at current power. The screen returned zero
FDR significant hits in the 100,204 case cohort. Both statements are true
simultaneously and neither should be reported without the other.

---

## 4. Result 2: triangulation across all three tests

Taxa whose effect points in the same direction in FinnGen colorectal,
Fernandez-Rozadilla colorectal, and FinnGen pancreatic cancer.

| Quantity | Value |
|---|---|
| Taxa present in all three tests | 210 |
| Consistent direction in all three | 82, or 39.0 percent |
| Naive expectation under independence | 52.5, or 25 percent |

**No p-value is reported for this excess, deliberately.** The three tests are
not independent. The two colorectal analyses share the same MiBioGen exposure
instruments, and taxa are structurally correlated through the taxonomic
hierarchy, which the WP0 nesting analysis showed can vary discovery counts up
to 5.7 fold. A binomial test against 25 percent would be anticonservative by
an unquantified amount. The excess is reported as a descriptive quantity only.

---

## 5. Result 3: Alistipes, corrected

**The previously recorded FinnGen colorectal value was wrong.** The record
carried OR 0.83 at p 0.026. No such value appears in any column of that row,
nor in any related Rikenellaceae lineage row.

| Test | Cases | OR | p (IVW) | p (weighted median) | FDR | Nominal |
|---|---|---|---|---|---|---|
| Colorectal, FinnGen R12 | 11,790 | **0.96** | **0.578** | 0.44 | 0.97 | **No** |
| Colorectal, Fernandez-Rozadilla | 100,204 | 0.93 | 0.048 | 0.278 | 0.67 | IVW only |
| Pancreatic, FinnGen | 3,139 | 0.68 | 0.017 | 0.035 | 0.89 | Yes, both |

Corrected statement, to be used in all downstream material:

> Genus *Alistipes* points in the protective direction in three of three
> independent tests and reaches nominal significance in two of three. The
> FinnGen colorectal test is null. The colorectal replication result is
> significant on inverse variance weighting only and does not hold on weighted
> median, with a divergent Egger estimate (OR 0.73, p 0.038, pleiotropy p
> 0.088). Only the pancreatic result holds on both IVW and weighted median
> with no pleiotropy flag. No test survives FDR correction in any cohort.
> This is a lead requiring species level instrumentation, not a finding.

The original honest framing, lead and not finding, is unaffected. Alistipes
remains the top ranked taxon on triangulation and the only one nominal in more
than one test. What changes is that the claim of three of three nominal
significance is withdrawn.

Biological rationale is unchanged: *Alistipes* is a butyrate adjacent
Bacteroidetes genus, and butyrate producers are the most repeatedly proposed
protective lineage in gastrointestinal cancer.

---

## 6. Figure

**Figure 7. The causal signal is real but diffuse.** (a) Log odds ratio in
FinnGen R12 against log odds ratio in the 100,204 case independent cohort,
210 gut taxa, sign recovered. Concordant taxa in blue, discordant in orange;
*Alistipes* circled. (b) Directional agreement under all three tie rules
against the 50 percent noise expectation, showing the excess does not depend
on the rule.

---

## 7. Outputs

| File | Content |
|---|---|
| `wp1_cross_consistency.py` | Canonical module, assertion checked, exits non-zero on drift |
| `cross_consistency.json` | Machine readable results including full tie sensitivity and top 10 triangulated taxa |
| `MR_cross_cohort_recon.tsv` | Per taxon reconstruction with b, se, stored p, reconstructed p, log10 delta, QC flag |
| `Figure7_cross_cohort.png` | 300 dpi, Okabe-Ito colorblind safe palette |

## 8. Residual limitations

1. The FinnGen colorectal arm is a **reconstruction**, not a re-derivation. It
   is validated against the stored p-value but the raw beta and standard error
   were never written. The definitive fix is to re-run the FinnGen colorectal
   MR with full precision output, which requires an OpenGWAS JWT and is not
   executable in this environment.
2. Phylum Cyanobacteria is excluded from the agreement statistic.
3. Triangulation excess is descriptive only, for the non independence reasons
   in section 4.
4. All limitations of the parent MR analyses carry through: liberal p < 1e-5
   instrument threshold, distance based rather than LD based clumping,
   palindromic SNPs dropped without frequency resolution, partial ancestry
   mismatch, no colocalization, no reverse direction MR.
