# MR Pipeline Stage 2 Validation

**Order:** RLO_MR_OUTCOME-STAGE2-VALIDATION_v1.0
**Program:** RO-2026-008, Causal Microbial Oncology
**Executed:** 2026-08-07
**Host:** escapepod-2.local, macOS arm64, Python 3.14, numpy 2.4.4, scipy 1.17.1
**Artefact under test:** `src/wp1_mr_pipeline.py` (committed at 2a3b4d6)

---

## 1. Plain answer

**The committed pipeline does not reproduce the committed result tables.**

Six distinct causes were isolated. Four are systematic and affect all 211 taxa in
all three outcomes. Two are conditional and affect a small number of taxa. A
seventh structural defect prevents the documented Stage 2 command from executing
at all.

Per RLO Section 8.7 this constitutes a PASS. The deliverable is a truthful
answer, not a confirmation.

---

## 2. Provenance of outcome data

| Outcome | Source | URL | Retrieved (UTC) | Bytes | sha256 |
|---|---|---|---|---|---|
| colorectal | FinnGen R12 `C3_COLORECTAL_EXALLC`, 11,790 cases / 378,749 controls | `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_COLORECTAL_EXALLC.gz` | 2026-08-07T14:07:35Z | 814,842,443 | `b42770f3470ec57cf5540b0e2d40f55f19f28865e04cdf6e501e30c499aafcac` |
| colorectal_100k | GWAS Catalog `GCST90129505`, Fernandez-Rozadilla 2023, 100,204 cases | `http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/harmonised/GCST90129505.h.tsv.gz` | 2026-08-07T14:26:35Z | 640,903,823 | `000b34b567c0540d13ce859476cfddcd9bafafac82d7eebc48a1bd3c74d8affc` |
| pancreatic | FinnGen R12 `C3_PANCREAS_EXALLC`, 3,139 cases / 378,749 controls | `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_PANCREAS_EXALLC.gz` | 2026-08-07T14:27:55Z | 810,871,731 | `e12eb6a465960bf7d61d0ca7e76729942b43077cfc17cc0708cad31d1df398ee` |
| exposure | MiBioGen `MBG.allHits.p1e4.txt`, 122,111 rows, 211 taxa | `https://molgenis26.gcc.rug.nl/downloads/MiBioGen/MBG.allHits.p1e4.txt` | 2026-08-07 | 18,329,457 | `37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e` |

The exposure sha256 matches the value recorded in the RLO exactly. RLO item A1
recorded a 404 for `finngen_R12_C3_COLORECTAL.gz`. The cause was the phenocode
and the path. The correct phenocode carries the `_EXALLC` suffix and the file
sits under `summary_stats/release/`, both confirmed against the live manifest
`finngen_R12_manifest.tsv` rather than assumed.

---

## 3. Stage 1 confirmation

Stage 1 was rerun on this host before any Stage 2 work.

```
taxa checked                 : 211
selected >= committed        : 211/211
selected == committed        : 52/211
median selected/committed    : 1.11
exit code                    : 0
```

The environment reproduces Stage 1 exactly as recorded in the RLO.

---

## 4. Structural defect found before any numerical comparison

`wp1_mr_pipeline.py` commits `load_exposure`, `select_instruments`, `ivw`,
`mr_egger`, `weighted_median`, `cochran_q` and `bh`. It does not commit a
working Stage 2. The `run` subcommand is a stub:

```python
print("Stage 2 is implemented but NOT validated end to end. ...", file=sys.stderr)
return 2
```

It never loads an outcome file, never harmonises, never loops over taxa and
never writes output. RLO item B1 specifies
`wp1_mr_pipeline.py run --exposure ... --outcome ... --out ...`. That command
cannot execute. Acceptance criterion 2 therefore cannot be satisfied by the
artefact as committed.

Validation proceeded by importing the committed functions unchanged into a
separate harness, `src/wp1_stage2_validate.py`, which supplies the
harmonisation and the per taxon loop that `run` never wired up. The committed
pipeline was not modified, in accordance with RLO Section 12.

---

## 5. Systematic discrepancies

Each was isolated by toggling it alone and measuring reproduction against the
committed tables across all 211 taxa.

| ID | Discrepancy | Committed pipeline | Table generating code | Consequence |
|---|---|---|---|---|
| D1 | Clumping window | 500 kb (`CLUMP_WINDOW`) | 1 Mb | different instrument sets for most taxa |
| D2 | Palindromic drop order | before clumping | after clumping, during harmonisation | changes which SNP wins each clump |
| D3 | MR-Egger orientation | none | exposure effects oriented positive | `p_egger` wrong for essentially every taxon |
| D4 | IVW standard error | fixed effect | random effects, scaled by sqrt(max(1, Q/df)) | `p_ivw` wrong wherever heterogeneity is present |

D3 is the most consequential for interpretation. Without orienting the exposure
effects positive, the MR-Egger slope and intercept are not the quantities the
method defines, and the reported pleiotropy test is not the Egger intercept
test. Reproduction of `p_egger` under committed parameters is 1 to 2 taxa out
of 211 across the three outcomes, which is consistent with chance.

D4 explains a specific signature seen during diagnosis: for taxa where the
instrument set happened to match, `b_ivw` and Cochran Q reproduced to the last
digit while `p_ivw` did not. That is a standard error difference, not an
estimator difference. The committed IVW point estimate and the committed
Cochran Q are correct.

---

## 6. Agreement summary, all 211 taxa, all three outcomes

Variant A applies the committed pipeline parameters exactly. Variant B applies
the faithful reconstruction with D1 through D4 corrected. Match tolerance is
1e-4 relative.

| Outcome | Variant | n_snp | p_ivw | p_egger | q_p | b (full precision) | se (full precision) |
|---|---|---|---|---|---|---|---|
| colorectal | A committed params | 70/211 | 32/211 | 1/211 | 70/211 | not comparable | not comparable |
| colorectal | B faithful | 208/211 | 207/211 | 207/211 | 208/211 | not comparable | not comparable |
| colorectal_100k | A committed params | 71/211 | 20/211 | 2/211 | 70/211 | not tested | not tested |
| colorectal_100k | B faithful | 211/211 | 209/211 | 209/211 | 209/211 | 209/211 | 209/211 |
| pancreatic | A committed params | 70/211 | 37/211 | 1/211 | 71/211 | not tested | not tested |
| pancreatic | B faithful | 208/211 | 208/211 | 208/211 | 208/211 | 208/211 | 208/211 |

The colorectal_100k row is the strongest evidence in this report. The faithful
reconstruction reproduces the committed full precision `b` and `se` for 209 of
211 taxa at float level. That establishes that the four systematic
discrepancies are a complete description of the difference between the
committed pipeline and the code that generated the tables, rather than a
partial one.

Per taxon comparison tables, one row per taxon for all 211 taxa, are written to:

| Outcome | Path |
|---|---|
| colorectal | `docs/findings/MR_pertaxon_comparison_colorectal.tsv` |
| colorectal_100k | `docs/findings/MR_pertaxon_comparison_colorectal_100k.tsv` |
| pancreatic | `docs/findings/MR_pertaxon_comparison_pancreatic.tsv` |

Each carries `committed_n`, `A_n`, `B_n`, `committed_p`, `A_p`, `B_p`,
`committed_b`, `B_b`, `committed_se`, `B_se`, `committed_p_egger`, `A_p_egger`,
`B_p_egger`, `committed_q_p`, `A_q_p`, `B_q_p`, `reproduces_A`, `reproduces_B`.

---

## 7. Conditional discrepancies, diagnosed individually

### D5. Multi allelic rsID row selection

FinnGen carries more than one row per rsID at multi allelic sites. Which row is
retained determines whether an instrument harmonises or is discarded as an
allele mismatch. Four such sites were identified in the colorectal outcome:

| rsID | Exposure alleles | Row 1 | Row 2 | Affected taxa |
|---|---|---|---|---|
| rs10805326 | G/A | G/A resolves | G/T mismatches | Peptostreptococcaceae, Intestinibacter |
| rs9580476 | C/T | T/C resolves | T/G mismatches | Porphyromonadaceae, Parabacteroides |
| rs2482038 | C/A | A/C resolves | A/T mismatches | Ruminiclostridium5 |
| rs2435610 | A/C | C/A resolves | C/T mismatches | Dialister |

Retaining the last row encountered leaves six taxa one instrument short in both
FinnGen outcomes. Retaining the first row recovers them. Reproduction improved
from 205/211 to 207/211 for colorectal and from 205/211 to 208/211 for
pancreatic on `p_ivw` under this change alone. The committed tables are
consistent with first row retention.

This is a genuine fragility, not a cosmetic one. The rule is undocumented in
the pipeline and the correct behaviour is to select the row whose allele pair
matches the exposure rather than to depend on file order.

### D6. Exposure rows carrying rsID "NA"

`MBG.allHits.p1e4.txt` contains rows whose `rsID` field is the literal string
`NA`. `select_instruments` applies no rsID validity filter, so such rows can be
selected as instruments. They can never be matched to an outcome by rsID and
are silently lost at harmonisation. Two taxa are affected in the colorectal
outcome, `order.Clostridiales.id.1863` and `class.Clostridia.id.1859`, each
carrying one `NA` instrument in its selected set.

The committed pipeline does not warn on this. An instrument that is selected,
counted toward instrument strength, and then silently dropped is a coverage gap
of the class this program treats as a defect rather than a nuisance.

### Residual non reproducing taxa, listed individually

| Outcome | Taxon | committed n | rebuilt n | committed p | rebuilt p | Relative difference | Diagnosis |
|---|---|---|---|---|---|---|---|
| colorectal | genus.Intestinibacter.id.11345 | 14 | 14 | 0.17930 | 0.17928 | 0.01% | numeric, just outside 1e-4 tolerance |
| colorectal | order.Clostridiales.id.1863 | 13 | 12 | 0.25737 | 0.28689 | 11.5% | one instrument lost, rsID `NA` (D6) plus mismatch at rs72915163 |
| colorectal | class.Clostridia.id.1859 | 12 | 11 | 0.49094 | 0.53542 | 9.1% | same as above |
| colorectal | genus..Eubacteriumfissicatenagroup.id.1437 | 10 | 9 | 0.68480 | 0.83032 | 21.3% | one instrument lost; taxon key differs between exposure and committed table |
| colorectal_100k | family.Oxalobacteraceae.id.2966 | 14 | 14 | 0.0015241 | 0.0015418 | 1.2% | counts match; one instrument differs in effect estimate |
| colorectal_100k | genus.Intestinibacter.id.11345 | 13 | 13 | 0.59953 | 0.60162 | 0.3% | as above |
| pancreatic | genus..Eubacteriumfissicatenagroup.id.1437 | 10 | 9 | 0.42438 | 0.83656 | 97.1% | one instrument lost, as colorectal |
| pancreatic | order.Clostridiales.id.1863 | 13 | 12 | 0.86188 | 0.91631 | 6.3% | rsID `NA` (D6) |
| pancreatic | class.Clostridia.id.1859 | 12 | 11 | 0.89858 | 0.95732 | 6.5% | rsID `NA` (D6) |

The same three taxa fail in both FinnGen outcomes, which confirms the cause is
deterministic and instrument side rather than outcome side.

---

## 8. Schema defects in the committed tables

The three committed tables use three different schemas.

| Table | Columns |
|---|---|
| `MR_results_colorectal.tsv` | `taxon, n_snp, OR (95% CI), p_ivw, fdr_ivw, p_egger, egger_intercept_p, p_wmedian, q_p` |
| `MR_results_colorectal_100k.tsv` | `taxon, n, OR, b, se, p, qp, or_eg, p_eg, pleio_p, p_wm, fdr` |
| `MR_results_pancreatic.tsv` | `taxon, n_snp, or_ivw, b_ivw, se_ivw, p_ivw, q, q_p, or_egger, p_egger, egger_intercept, egger_intercept_p, or_wmedian, p_wmedian, fdr_ivw` |

`MR_results_colorectal.tsv` stores the effect only as the rounded display string
`0.79 (0.71-0.89)`. It carries no `b` and no `se`. Full precision reproduction
of that table is therefore not checkable, which is the residual of CORRECTIONS
C7 that RLO item B3 targets. The other two tables do carry full precision `b`
and `se` and were checked at float level.

RLO item B3 is satisfied in the following sense. The rebuilt colorectal table
now exists with full precision `b` and `se` at
`results/MR_rebuilt_colorectal.tsv`, produced under the faithful reconstruction
that reproduces 207 of 211 taxa. It is a reconstruction, not a recovery of the
original floats, and it is labelled as such. The original floats do not exist
anywhere in the repository and cannot be recovered from a rounded string.

---

## 9. What this does and does not license

| Statement | Status |
|---|---|
| The committed MR result tables are wrong | NOT established. This order tested the pipeline against the tables, not the tables against the data |
| The committed pipeline reproduces the tables | Refuted. Four systematic and two conditional discrepancies |
| The method actually used to generate the tables is now identified | Established, to 209/211 at float level on the one table carrying full precision values |
| Any published claim must change | NOT established by this order. RLO Section 12 forbids claim edits here |

The distinction matters. The faithful reconstruction reproduces the committed
numbers, which means the tables are internally consistent with a coherent and
conventional MR implementation using 1 Mb clumping, palindromic removal at
harmonisation, random effects IVW and oriented MR-Egger. The defect is that the
committed pipeline is not that implementation. The derivation path was missing
and, where it existed, it disagreed with the artefact it claimed to document.

## 10. Artefacts produced

| Artefact | Path |
|---|---|
| Rebuilt table, colorectal, faithful reconstruction | `results/MR_rebuilt_colorectal.tsv` |
| Rebuilt table, colorectal 100k, faithful reconstruction | `results/MR_rebuilt_colorectal_100k.tsv` |
| Rebuilt table, pancreatic, faithful reconstruction | `results/MR_rebuilt_pancreatic.tsv` |
| Rebuilt table, colorectal, committed pipeline parameters | `results/MR_rebuilt_colorectal_committedparams.tsv` |
| Rebuilt table, colorectal 100k, committed pipeline parameters | `results/MR_rebuilt_colorectal_100k_committedparams.tsv` |
| Rebuilt table, pancreatic, committed pipeline parameters | `results/MR_rebuilt_pancreatic_committedparams.tsv` |
| Per taxon comparisons | `docs/findings/_pertaxon_{colorectal,colorectal_100k,pancreatic}.tsv` |
| Stage 2 harness | `src/wp1_stage2_validate.py` |
| Comparison and diagnosis | `src/wp1_stage2_compare.py`, `src/wp1_stage2_diagnose.py` |
| Outcome downloads | `outcome/` (gitignored, not committed) |

Both variants are emitted so the artefact under test and the reconstruction that
reproduces the tables can be compared directly. `--faithful` selects the
reconstruction; the default reproduces the committed pipeline exactly.

Verification of the faithful tables against the committed tables:

| Outcome | n_snp | p_ivw | b | se |
|---|---|---|---|---|
| colorectal | 208/211 | 207/211 | not in committed table (C20) | not in committed table (C20) |
| colorectal_100k | 211/211 | 209/211 | 209/211 | 209/211 |
| pancreatic | 208/211 | 208/211 | 208/211 | 208/211 |

All rebuilt tables carry full precision `b` and `se` as floats per RLO
Section 5. No display strings are used to record an effect.

## 11. Acceptance criteria

| # | Criterion | Status |
|---|---|---|
| 1 | URL, date, size, sha256 for every outcome file | Met, Section 2 |
| 2 | `wp1_mr_pipeline.py run` completes on all three outcomes | NOT MET. The subcommand is a stub that exits 2. Documented as the finding, Section 4 |
| 3 | Rebuilt tables carry full precision b and se | Met, Section 10 |
| 4 | Per taxon comparison for all 211 taxa | Met, Section 6 |
| 5 | Every discrepancy individually listed and diagnosed | Met, Sections 5 and 7 |
| 6 | Report states plainly whether the tables reproduce | Met, Section 1 |
| 7 | Non reproduction is a PASS | Invoked |
| 8 | No repo claim changed on unvalidated output | Met. No README, manuscript, poster or figure was modified |
