# RO-2026-008 | WP0 | SIA_001
## Analysis 7: taxonomic nesting redundancy (real TCMA data)

Taxonomic hierarchies produce structurally collinear features: a genus column
is by construction a function of its species columns. Two consequences follow
that are not consistently controlled in the tumour microbiome literature.

### A. Discovery counts depend on the redundancy rule, and vary up to 5.7-fold

| Analysis | Reported | r>0.99 | r>0.95 | r>0.90 | r>0.80 | Ancestor rule |
|---|---|---|---|---|---|---|
| STAD, paired | 8 | 3 | 3 | 2 | 2 | 2 |
| HNSC, unpaired | 51 | 44 | 43 | 42 | 31 | 9 |

Same data, same p-values, same FDR threshold. The number of "findings" ranges
from 2 to 8 in STAD and from 9 to 51 in HNSC purely as a function of how
nested ranks are collapsed.

Neither extreme is correct. The pure ancestor rule over-collapses: it merges
41 distinct organisms into a single "Actinobacteria" clade, discarding real
genus-level distinctions between Rothia, Actinomyces and Veillonella. The
r>0.99 rule under-collapses: it retains near-duplicate strain entries
(Rothia sp. HMSC071B01, HMSC078H08, HMSC072B03 and so on) as separate
discoveries. The defensible range sits between.

### B. A substantial fraction of the feature space is structurally redundant

Proportion of taxa surviving as independent features (first 500 taxa per cohort):

| Cancer | Nominal | r>0.99 | r>0.95 | r>0.90 | Redundant at r>0.90 |
|---|---|---|---|---|---|
| COAD | 500 | 418 | 388 | 368 | 26% |
| STAD | 452 | 359 | 329 | 306 | 32% |
| HNSC | 500 | 402 | 384 | 366 | 27% |
| ESCA | 431 | 310 | 244 | 202 | 53% |

Between roughly one quarter and one half of the feature space carries no
information beyond an ancestor rank.

### C. Implication for multiple-testing control

BH-FDR assumes independence or positive regression dependence. Perfectly
nested features do not violate PRDS in a way that invalidates BH, so the
error rate is not necessarily wrong. The problem is interpretive rather than
inferential: the nominal test count overstates the number of independent
biological hypotheses, and the discovery count overstates the number of
independent findings. Reported taxa counts are therefore not comparable
across studies unless the collapsing rule is stated, and it rarely is.

### Recommendation for the field
1. State the redundancy rule explicitly in Methods.
2. Report discoveries as independent clades, not raw taxon counts.
3. Report the effective feature count alongside the nominal count.

### Honest limitation
The claim that this rule is "rarely stated" is our impression from the
anchor corpus, not a systematic survey. A structured survey of published
tumour microbiome papers would be required to substantiate it, and is not
yet done.
