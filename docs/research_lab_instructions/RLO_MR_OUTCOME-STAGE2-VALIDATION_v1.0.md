# RLO_MR_OUTCOME-STAGE2-VALIDATION_v1.0

**Issued:** 2026-08-07
**Program:** RO-2026-008, Causal Microbial Oncology
**Repo:** julian-borges-md/tumor-microbiome-causality-audit @ 2a3b4d6
**Blocker status:** L2 partially closed. Stage 1 validated. Stage 2 blocked on data acquisition only.

---

## 1. What I Need

| # | Output | Path |
|---|---|---|
| 1 | Outcome summary statistics, harmonised to the committed instrument set | `outcome/` (gitignored, large) |
| 2 | Validation report proving `wp1_mr_pipeline.py` Stage 2 reproduces the committed tables, or documenting exactly where it does not | `docs/findings/MR_PIPELINE_VALIDATION.md` |

---

## 2. Why It Exists

`src/wp1_mr_pipeline.py` commits the full MR method. Stage 1, exposure instrument selection, is validated against all 211 committed `n_snp` counts with zero deficits. Stage 2, outcome harmonisation plus IVW, MR-Egger, weighted median, Cochran Q and BH-FDR, is written but has never been run against real outcome data. Its output must not be used for any claim.

This is the last structural gap in the program. Until it closes, the causal arm has results without a derivation path, which is the exact defect class that produced CORRECTIONS C5.

---

## 3. What Has Been Done

| Item | State |
|---|---|
| Exposure data | MiBioGen `MBG.allHits.p1e4.txt`, 122,111 rows, 211 taxa, sha256 `37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e`. Public, no credentials |
| Instrument parameters | **Recovered, not assumed.** p<1e-5, 500kb distance clumping, palindromic dropped. Zero deficits over 211 taxa, median selected/committed 1.11. See CORRECTIONS C12 |
| Estimators | Implemented in `wp1_mr_pipeline.py` |
| Committed truth tables | `results/MR_results_colorectal.tsv`, `MR_results_colorectal_100k.tsv`, `MR_results_pancreatic.tsv` |
| Outcome data | **Not obtained.** This is the entire blocker |

---

## 4. Remaining Work

| ID | Item | Notes |
|---|---|---|
| A1 | FinnGen R12 colorectal cancer summary statistics | Phenocode unknown. `finngen_R12_C3_COLORECTAL.gz` returned 404 on the public bucket. Locate the R12 manifest; do not guess |
| A2 | GCST90129505, Fernandez-Rozadilla 2023, 100,204 cases | GWAS Catalog FTP, harmonised release preferred |
| A3 | FinnGen R12 pancreatic cancer summary statistics | Same manifest as A1 |
| A4 | Fallback: register an OpenGWAS JWT | api.opengwas.io returns 401 unauthenticated. Valid 14 days. Strip the `1-` UI prefix before use |
| B1 | Run Stage 2 per outcome | `wp1_mr_pipeline.py run --exposure ... --outcome ... --out ...` |
| B2 | Per-taxon comparison, rebuilt versus committed | All 211 taxa, not sampled |
| B3 | Resolve L3: re-export FinnGen colorectal with full-precision `b` and `se` | The committed table stores a rounded display string only |

---

## 5. Schema Per Entry

`results/MR_rebuilt_{outcome}.tsv`, one row per taxon:

```yaml
taxon:             str      # e.g. genus.Alistipes.id.968
n_snp:             int      # post-harmonisation instrument count
b_ivw:             float    # log OR, FULL PRECISION
se_ivw:            float
p_ivw:             float
or_ivw:            float
lci:               float
uci:               float
fdr_ivw:           float    # BH across all taxa in this outcome
b_egger:           float
se_egger:          float
p_egger:           float
egger_intercept_p: float
b_wmedian:         float
se_wmedian:        float
p_wmedian:         float
q_stat:            float
q_p:               float
```

**Non-negotiable:** `b` and `se` are stored as floats. Never write a formatted `OR (LCI-UCI)` display string as the only record of an effect. That is CORRECTIONS C7 and it cost a full reconstruction.

---

## 6. Validation Report Requirements

1. Provenance: URL, download date, file size, sha256 for every outcome file
2. Per-taxon comparison table, rebuilt versus committed, all 211 taxa
3. Agreement summary: how many reproduce within tolerance, how many do not
4. **Every discrepancy listed individually with a diagnosis.** Not an aggregate
5. One plain sentence stating whether the committed tables reproduce

---

## 7. Databases to Query

| Source | What | Access notes |
|---|---|---|
| FinnGen R12 | Colorectal, pancreatic | Public. Locate the release manifest first |
| GWAS Catalog | GCST90129505 | Public FTP, harmonised release |
| OpenGWAS REST | Fallback | JWT required. GET `/api/gwasinfo?id={id}` is a query param, not a path. POST `/api/associations` uses field `variant`, not `rsid` |
| MiBioGen | Already obtained | No action |

---

## 8. Acceptance Criteria

1. Every outcome file has URL, date, size and sha256 recorded
2. `wp1_mr_pipeline.py run` completes without error on all three outcomes
3. Rebuilt tables carry full-precision `b` and `se`, no display strings
4. Per-taxon comparison reported for all 211 taxa
5. Every discrepancy individually listed and diagnosed
6. The report states plainly whether the committed tables reproduce
7. **If they do not reproduce, that is a PASS for this order.** The deliverable is a truthful answer, not a confirmation. A discrepancy found here is worth more than one found by a reviewer
8. No repo claim is changed on the basis of unvalidated Stage 2 output

---

## 9. Output File Locations

| Artefact | Path |
|---|---|
| Rebuilt tables | `results/MR_rebuilt_{colorectal,colorectal_100k,pancreatic}.tsv` |
| Validation report | `docs/findings/MR_PIPELINE_VALIDATION.md` |
| Provenance | append to `docs/DATA_CHECKSUMS.txt` |
| Corrections | append to `docs/CORRECTIONS.md` if anything fails to reproduce |
| Raw downloads | `outcome/` — **gitignore, do not commit summary statistics** |

---

## 10. What This Unblocks

| Blocked item | Unblocked by |
|---|---|
| L2, MR pipeline uncommitted | B1 + B2 |
| L3, FinnGen reconstruction | B3 |
| Reconstruction step in `wp1_cross_consistency.py` | B3. The log-CI midpoint recovery becomes unnecessary |
| Second manuscript, the MR paper | All of the above. Not submittable while its pipeline is unvalidated |
| CORRECTIONS C7 residual | B3 |

---

## 11. Step Sequence

1. Locate the FinnGen R12 release manifest. Record exact phenocodes for colorectal and pancreatic cancer. Do not guess.
2. Download A1, A2, A3. Record sha256 before use.
3. If public routes fail, register an OpenGWAS JWT and use the REST fallback.
4. Run Stage 1 to confirm the environment reproduces: 211/211, zero deficits, exit 0.
5. Run Stage 2 per outcome.
6. Build the per-taxon comparison against the committed tables.
7. Write the validation report. List every discrepancy.
8. Append provenance to `DATA_CHECKSUMS.txt`.
9. If anything fails to reproduce, append a numbered correction to `CORRECTIONS.md`. Do not silently adjust the pipeline until it matches.
10. Commit. Do not push changes to any claim in README or on the poster without review.

---

## 12. Explicit Prohibitions

| Do not | Reason |
|---|---|
| Tune Stage 2 parameters until output matches the committed tables | That converts a validation into a curve fit. Report the mismatch |
| Commit downloaded summary statistics | Size, and redistribution terms vary by source |
| Update any claim in README, manuscript or poster | Requires review. This order produces evidence, not edits |
| Report aggregate agreement without the per-taxon table | CORRECTIONS C6 exists because a summary statistic hid a rule dependency |
| Treat a null or a mismatch as failure | A mismatch found here is the most valuable output this order can produce |
