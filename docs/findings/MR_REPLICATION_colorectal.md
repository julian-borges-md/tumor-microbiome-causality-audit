# RO-2026-008 | WP1 | Colorectal MR: replication in a large independent cohort

## Design
The FinnGen screen (11,790 cases) produced suggestive leads. To test whether they
were real, the same instruments were run against an INDEPENDENT, larger cohort:
Fernandez-Rozadilla et al. 2023 (Nat Genet), 100,204 cases / 154,587 controls,
GWAS Catalog GCST90129505. This is ~8.5x the FinnGen case count.

## Result: the leads did NOT replicate

| Taxon | FinnGen OR (p) | Fernandez-Rozadilla OR (p) | Replicated? |
|---|---|---|---|
| genus Sutterella | 0.81 (2.3e-3) | 1.02 (0.47) | NO |
| family Enterobacteriaceae | 1.46 (6.8e-3) | 1.01 (0.87) | NO |
| order Enterobacteriales | 1.46 (6.8e-3) | 1.01 (0.87) | NO |
| phylum Cyanobacteria (FinnGen FDR hit) | 0.79 (1.2e-4) | 1.10 (0.03), direction FLIPPED | NO, confirmed artifact |
| genus Bifidobacterium | 1.21 (2.1e-2) | 1.06 (0.15) | NO (was LCT-driven) |

## In the large cohort, no taxon survives correction
211 taxa tested against 100,204 cases. FDR-significant hits: ZERO. Strongest new
signals (family Oxalobacteraceae OR 0.93 p 1.5e-3; order Clostridiales OR 1.15
p 3.3e-3) sit at FDR 0.15, suggestive only, and given the fate of the FinnGen leads
should be treated as probable false positives pending their own replication.

## Interpretation
1. The FinnGen leads were false positives. Independent replication at higher power
   eliminated them. This is the intended function of replication and is the most
   trustworthy result in the program.
2. At the largest available colorectal cancer GWAS, gut bacteria show no robust
   causal effect on colorectal cancer risk surviving multiple-testing correction.
3. The OUTCOME is no longer the limiting factor: 100k cases still yield null. The
   bottleneck is the EXPOSURE. MiBioGen is a 2021 genus-level 16S GWAS (n=18,340)
   with weak instruments that cannot resolve species. The Fusobacterium nucleatum
   to colorectal cancer link operates at the species level and may be uninstrumentable
   with genus-level 16S data.

## The genuine next lever: species-level metagenomic microbiome GWAS
Species-level shotgun-metagenomic mGWAS now exist and could instrument specific
cancer-associated organisms:
- Tomofuji et al. 2023 (Cell Rep), Japanese, 423 microbial species, public resource.
- Swedish shotgun-metagenome GWAS, n=16,017 + Norwegian HUNT replication (n=12,652),
  species-level.
- Zhernakova et al. 2024 (Nature), Dutch, structural variants across 49 species.
Constraints: (a) ancestry matching (Japanese exposure needs East Asian cancer outcome);
(b) few instruments per species (microbiome is weakly heritable), so tests remain fragile;
(c) uncertain whether Fusobacterium nucleatum specifically has genome-wide instruments.
This is the only move that could genuinely change the current null.
