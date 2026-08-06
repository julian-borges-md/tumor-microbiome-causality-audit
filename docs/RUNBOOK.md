# RO-2026-008 | WP0 | Reproducibility Runbook
Version wp0_core-1.0.0

## 1. Environment
See ENVIRONMENT.txt. Verified on python 3.12.3, numpy 2.4.4, pandas 3.0.2,
scikit-learn 1.8.0, scipy 1.17.1, Linux x86_64.

## 2. Data acquisition
Source: The Cancer Microbiome Atlas, Dohlman et al., Cell Host & Microbe 2021
DOI: 10.7924/r4bk1j35s  (resolves to research.repository.duke.edu/record/126)

    mkdir -p tcma && cd tcma
    for f in README.txt metadata.zip WGS.zip; do
      curl -sL "https://research.repository.duke.edu/record/126/files/$f?version=1" -o "$f"
    done
    unzip -o -q metadata.zip -d meta && unzip -o -q WGS.zip -d wgs

Verify against DATA_CHECKSUMS.txt before proceeding:

    sha256sum -c DATA_CHECKSUMS.txt

## 3. Synthetic analyses (deterministic given config + seed)

    python3 wp0_core.py audit --seeds 0,1,2,3,4,5,6,7,8,9 --out audit.json
    python3 wp0_core.py sweep --seeds 0,1,2                --out sweep.json
    python3 wp0_core.py floor --seeds 0,1,2,3,4,5,6,7,8,9  --out floor.json

## 4. Real-data analyses with verification

    python3 wp0_tcma.py

Exits non-zero and prints DRIFT lines if any reported value fails to
reproduce within tolerance. Exit 0 means all values verified.

## 5. Tolerances
Accuracies: absolute 0.02. H. pylori paired difference: absolute 0.01.
log10 p-value: absolute 0.15. Sample and pair counts: exact.

## 6. Determinism notes
- Seeds are explicit arguments, never global state. Internal RNG is
  numpy default_rng(1_000_000 + seed).
- Classifier random_state and CV random_state are fixed at 0, so
  classifier variance is not a source of run-to-run drift; simulation
  seed is the only stochastic input.
- Consequence: dispersion reported across seeds reflects SIMULATION
  variance only, not classifier variance. This understates total
  uncertainty and is stated as a limitation.
- Real-data analyses are fully deterministic; no seed dependence.
