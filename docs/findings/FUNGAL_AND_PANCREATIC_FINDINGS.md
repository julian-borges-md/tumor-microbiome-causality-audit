# RO-2026-008 | WP2/WP1 | Fungal causation feasibility + bacterial MR to pancreatic cancer

## Part 1: Can Malassezia (and other fungi) be tested causally? Feasibility determination.

Motivation: Aykut et al. (Nature 2019) showed Malassezia drives pancreatic ductal
adenocarcinoma in mice via mannose-binding lectin (MBL) and complement C3 activation.
The causal question in humans requires genetic instruments for fungal abundance.

Finding: DIRECT MR OF MALASSEZIA IS NOT CURRENTLY FEASIBLE with accessible data.

| Resource | Malassezia included? | Accessible summary stats? | Power |
|---|---|---|---|
| PLOS Biology 2025 (first mycobiome GWAS, pbio.3003339) | NO (9 taxa: Aspergillus, Candida, Kazachstania, Saccharomycetaceae, Aspergillaceae, Pleosporales, Capnodiales, Saccharomycetales, candidate Saccharomycetales family) | Yes (open) | Very low: HMP n=125, many hits sub-genome-wide |
| Vita 2026 (65-taxon mycobiome GWAS) | YES (Malassezia colocalizes with TNK2-AS1 eQTL) | No accessible summary statistics located | Larger but unavailable |
| CARD9 locus | Determines Malassezia gut colonization | Effect size on Malassezia requires a mycobiome GWAS (unavailable) | n/a |

This is the fungal analog of the CYP19A1 constraint (zero brain cis-eQTLs): the direct
causal mechanism cannot be tested with data that exists in downloadable form. Malassezia
has, at best, one to two genome-wide loci in the single GWAS that includes it, and that
GWAS is not accessible. A direct test would be a fragile single-instrument Wald ratio with
no ability to test the exclusion restriction.

### The two real paths to test the fungal-pancreatic hypothesis
1. MECHANISM TEST: instrument MBL / mannose-binding lectin (MBL2 has large, well-characterized
   cis effects on protein level) and complement C3, and run MR to pancreatic cancer. Tests
   whether the lectin-complement arm that Malassezia exploits causally affects PDAC in humans.
   Requires plasma-proteome pQTL instruments (deCODE / UKB-PPP), a further data step.
2. WAIT FOR DATA: the 2026 mycobiome GWAS with Malassezia signal, once summary statistics
   become accessible, would permit a direct (still instrument-limited) test.

## Part 2: Executable causal test. Bacterial gut taxa to pancreatic cancer.

Since the pancreatic outcome and the bacterial pipeline were in hand, the causal test that
IS executable was run: 211 MiBioGen gut taxa against FinnGen pancreatic cancer
(C3_PANCREAS_EXALLC, 3,139 cases). This does not test fungi. It tests the same underlying
question for bacteria.

| Taxon | N SNP | OR (95% CI) | p (IVW) | FDR | Egger intercept p | Q p | LOO robust | Verdict |
|---|---|---|---|---|---|---|---|---|
| genus (unclassified) id.959 | 13 | 1.32 (1.14-1.52) | 1.9e-4 | 0.04 | 0.47 | 0.67 | Yes | FDR-significant, statistically CLEAN, but UNIDENTIFIED organism |
| genus LachnospiraceaeUCG004 | 14 | 0.68 (0.53-0.88) | 2.9e-3 | 0.31 | 0.55 | 0.47 | Yes | Suggestive protective; butyrate producer; biologically coherent |
| genus Alistipes | 12 | 0.68 (0.50-0.93) | 1.7e-2 | 0.89 | 0.94 | 0.97 | No | Fragile |
| genus Butyricicoccus | 8 | 1.45 (1.07-1.97) | 1.7e-2 | 0.89 | 0.54 | 0.65 | No | Fragile, single-SNP driven |

### Interpretation
The single FDR-surviving hit is statistically the cleanest in either cancer screen: no
pleiotropy flag, no heterogeneity, concordant weighted median, robust to leave-one-out, and
its instruments show no obvious pleiotropy through a known pancreatic pathway. But it is an
UNCLASSIFIED MiBioGen genus, so the organism cannot be named, which severely limits
interpretation and any translational use. The one biologically coherent signal
(LachnospiraceaeUCG004, a butyrate producer, protective) is suggestive only.

### Honest bottom line
For bacteria, one clean but unidentifiable genus causally raises pancreatic cancer risk at
this power, and one plausible butyrate-producer lineage is suggestively protective. For
fungi, and specifically Malassezia, the direct causal test cannot yet be run. The
mechanism test (MBL and complement to pancreatic cancer) is the tractable way to probe the
Aykut pathway and is the recommended next step.

## Shared limitations with the colorectal analysis
Instrument threshold 1e-5, distance-based clumping, palindromic SNPs dropped, single Finnish
outcome cohort, no colocalization yet, no reverse MR yet. See MR_FINDINGS_colorectal.md.
