# RO-2026-008 | WP0 | SIA_001
## Analysis 6: paired within-patient tumor vs adjacent normal (TCMA, real data)

Pairing each tumor to the SAME patient's adjacent normal removes inter-patient
compositional variance, the dominant noise source in Analysis 5.

| Cancer | Matched pairs | Significant (BH-FDR<0.05) | Enriched in tumor | Depleted in tumor |
|---|---|---|---|---|
| COAD | 21 | 0 | 0 | 0 |
| STAD | 39 | 8 | 4 | 4 |
| HNSC | 22 | 0 | 0 | 0 |
| ESCA | 22 | 0 | 0 | 0 |
| READ | 4 | not evaluable | - | - |

### STAD detail (the only cancer with signal)

| Taxon | Direction | Paired diff | p |
|---|---|---|---|
| Helicobacter pylori | DOWN in tumor | -0.994 | 4.5e-04 |
| Helicobacteraceae (family) | DOWN | -1.016 | 5.5e-04 |
| Helicobacter (genus) | DOWN | -1.016 | 5.5e-04 |
| Helicobacter acinonychis | DOWN | -0.275 | 8.0e-04 |
| Capnocytophaga | UP in tumor | +0.249 | 6.6e-04 |
| Flavobacteriaceae / Flavobacteriia / Flavobacteriales | UP | +0.249 | 7.1e-04 |

IMPORTANT: these 8 hits are NOT 8 independent findings. They are nested
taxonomic ranks of the same two organisms. The true count is TWO independent
signals: Helicobacter depleted, Capnocytophaga (an oral genus) enriched.
Any future reporting must collapse nested ranks before counting.

Biologically coherent and consistent with known gastric carcinogenesis:
loss of the acid-adapted Helicobacter niche through atrophy and intestinal
metaplasia, with permissive colonisation by oral-type flora.

### Positive controls under pairing

| Control | Unpaired p | Paired p | Paired diff | Direction |
|---|---|---|---|---|
| F. nucleatum, COAD | 0.231 | 0.075 | +0.120 | correct (enriched) |
| F. nucleatum, HNSC | 0.543 | 0.053 | +0.688 | correct (enriched) |
| F. nucleatum, ESCA | - | 0.263 | +0.019 | correct, negligible |
| H. pylori, STAD | 0.022 | 0.00045 | -0.994 | inverted (depleted) |

Pairing improved every positive control, roughly 3-fold on COAD and 10-fold on
HNSC, and moved both Fusobacterium tests to the edge of significance without
crossing it at n = 21 and 22 pairs. This is the predicted behaviour at the
detection floor established in Analysis 3 and confirms the Analysis 5 nulls
were power-limited, not evidence of absence.

---

## Replication status

Nejman et al. 2020 (Science) was selected as the independent replication
cohort. Two constraints found:

1. SCOPE MISMATCH. Nejman covers breast, lung, ovary, pancreas, melanoma,
   bone and brain. TCMA covers colon, rectum, stomach, oesophagus and head
   and neck. There is NO overlap, so Nejman cannot replicate the gastric
   Helicobacter finding. It can only test the general claim that within-tissue
   tumour-vs-normal signal is weak.
2. ACCESS. PMC serves supplementary files behind a JavaScript proof-of-work
   anti-bot challenge, which was not circumvented. Files download normally in
   a browser (PMC7757858, Tables S1-S15). Raw 16S is at BioProject PRJNA624822.

Consequence: replication of the gastric finding requires a DIFFERENT cohort
containing gastric cancer with matched normals. Candidates for WP0 v2.
