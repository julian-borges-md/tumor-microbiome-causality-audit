# RO-2026-008 | WP1 | Two-sample MR: gut microbiome to colorectal cancer

## Design
- Exposure: MiBioGen gut microbiome GWAS, 211 taxa, 18,340 individuals (Kurilshikov 2021).
- Outcome: FinnGen R12 colorectal cancer, C3_COLORECTAL_EXALLC, 11,790 cases / ~300k controls.
- Instruments: p < 1e-5, distance-clumped to 1 Mb independence. 2,922 instruments, all 211 taxa with >=3, median F 21.
- Harmonization: matched on rsID, 448 palindromic SNPs dropped (no exposure EAF to resolve strand), 10 allele mismatches dropped, 2,370 usable.
- Estimators: IVW (random effects) primary; MR-Egger with intercept pleiotropy test, weighted median, Cochran Q. BH-FDR across 211 taxa.
- LD clumping was distance-based (no LD reference panel available); this is an approximation to r2-based clumping and is a stated limitation.

## Headline result
Across 211 gut taxa, ONE survives FDR correction, and it is most likely a pleiotropic artifact. No taxon shows a robust, biologically coherent, FDR-surviving causal effect on colorectal cancer risk at this power.

| Taxon | N SNP | OR (95% CI) | p (IVW) | FDR | Egger intercept p | Q p | Verdict |
|---|---|---|---|---|---|---|---|
| phylum Cyanobacteria | 8 | 0.79 (0.71-0.89) | 1.2e-4 | 0.02 | 0.05 | 0.68 | FDR-significant but LIKELY ARTIFACT |
| genus Sutterella | 12 | 0.81 (0.71-0.93) | 2.3e-3 | 0.24 | 0.49 | 0.45 | Cleanest signal; suggestive only |
| family Enterobacteriaceae | 7 | 1.46 (1.11-1.93) | 6.8e-3 | 0.25 | 0.36 | 0.03 | Biologically coherent; suggestive; mild heterogeneity |
| genus Bifidobacterium | 15 | 1.21 (1.03-1.41) | 2.1e-2 | 0.44 | 0.37 | 0.01 | ARTIFACT: driven by LCT locus (see below) |

## Why the FDR hit is likely not causal
The Cyanobacteria signal is statistically robust in the narrow sense: 7 of 8 per-SNP Wald ratios point protective, it survives leave-one-out (max p 1.4e-3), and estimators are concordant. But:
1. Its Egger intercept sits at p = 0.05, flagging directional pleiotropy.
2. One instrument (rs2314810) is at NOS2, inducible nitric oxide synthase, a gene directly implicated in colorectal carcinogenesis. An instrument that affects the outcome through a known cancer pathway violates the exclusion restriction.
3. "Gut Cyanobacteria" in 16S data commonly reflects dietary chloroplast DNA or the non-photosynthetic class Melainabacteria, so the exposure itself is of doubtful biological identity.
Taken together, this is most parsimoniously read as pleiotropy, not a causal effect of a bacterium.

## Pipeline validation: the LCT catch
The Bifidobacterium risk signal (OR 1.21) collapses entirely when a single SNP, rs182549, is removed (OR falls to 1.08, p rises to 0.27). rs182549 is at the LCT lactase locus, the canonical Bifidobacterium instrument, which is confounded with dairy intake and is heavily pleiotropic. The analysis correctly identifies this known trap through leave-one-out, which is evidence the rigor is working rather than rubber-stamping.

## The two signals worth following (suggestive, not confirmed)
- genus Sutterella, protective, OR 0.81. The cleanest result in the screen: no pleiotropy flag, no heterogeneity, LOO-robust, instruments show no obvious pleiotropy (autophagy and signaling genes, no diet loci). Does not survive FDR.
- family Enterobacteriaceae, risk-increasing, OR 1.46. Biologically coherent, since Enterobacteriaceae includes colorectal-cancer-associated pathobionts (e.g., certain E. coli). LOO-robust but mild heterogeneity, and does not survive FDR.

## Honest bottom line
This is a causal analysis, with direction fixed by germline genetics, so even the predominantly null result is a causal statement: at 11,790 cases, most gut taxa do not measurably cause colorectal cancer, the one taxon passing correction is best explained by pleiotropy, and the biologically plausible candidates are suggestive but unconfirmed. This is more informative than any cross-sectional abundance analysis, and it is not the discovery of a causal driver.

## Limitations
1. Instrument threshold p<1e-5 is liberal (standard in microbiome MR) and admits weaker instruments; F is bounded below by the threshold, so the reported F does not independently establish strength.
2. Distance-based clumping approximates LD-based clumping; residual correlation between instruments is possible.
3. Palindromic SNPs dropped without frequency resolution, reducing instruments.
4. Single outcome cohort (FinnGen, Finnish ancestry); MiBioGen is multi-ancestry but European-majority. Ancestry mismatch is partial.
5. No colocalization performed yet; a shared causal variant has not been confirmed for any hit.
6. Reverse-direction MR not yet run.
7. One cancer only. Other sites (gastric, given the H. pylori result) not yet tested.
