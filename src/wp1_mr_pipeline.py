#!/usr/bin/env python3
"""
RO-2026-008 | WP1 | Two-sample Mendelian randomization pipeline
wp1_mr_pipeline-1.0.0

Closes limitation L2 in part. The MR arm previously existed only as result
tables and plots: instrument selection, harmonisation and every estimator were
run interactively and never committed. That is the same defect class that
produced the Alistipes error (CORRECTIONS C5), one level up.

Scope of THIS module: Stage 1 only.

  STAGE 1  EXPOSURE. Instrument selection from MiBioGen, publicly
           downloadable, no credentials. Parameters were RECOVERED from the
           committed result tables rather than assumed; see
           validate_instruments() and docs/CORRECTIONS.md C12.

  STAGE 2  OUTCOME. NOT in this module. An earlier version carried a `run`
           subcommand that printed a notice and returned 2: no outcome
           loading, no harmonisation, no taxon loop, no writer. It was
           documented as "implemented", which was false, and it was removed
           (docs/CORRECTIONS.md C13, C21).

           The canonical Stage 2 implementation is src/wp1_stage2_validate.py.

Estimators (ivw, mr_egger, weighted_median, cochran_q, bh) live here and are
imported by the Stage 2 harness. Note that as written they carry four
diagnosed defects relative to the code that produced the published tables:
fixed-effect rather than random-effects IVW standard errors, MR-Egger without
exposure orientation, a 500kb rather than 1Mb clumping window, and palindromic
removal before rather than during harmonisation. See docs/CORRECTIONS.md
C14-C17. A flag audit found these change 0/211 published Figure 6
classifications, so they are code defects, not result errors.

Usage
-----
    # Stage 1, runs today, no credentials
    python3 src/wp1_mr_pipeline.py validate-instruments \\
        --exposure /path/to/MBG.allHits.p1e4.txt \\
        --committed results/MR_results_colorectal.tsv

    # Stage 2 lives elsewhere
    python3 src/wp1_stage2_validate.py --exposure ... --outcome-hits ... \\
        --fmt finngen --out results/MR_rebuilt_colorectal.tsv --faithful

Exit codes
----------
    0  validation passed
    1  DRIFT
    2  missing input
"""

import argparse
import csv
import math
import os
import statistics as st
import sys

import numpy as np
from scipy import stats

# Recovered instrument-selection parameters. See validate_instruments().
P_THRESHOLD = 1e-5
CLUMP_WINDOW = 500_000          # base pairs, distance-based
DROP_PALINDROMIC = True
PALINDROMIC = {frozenset("AT"), frozenset("CG")}

# Assertions for stage 1, over the 211 taxa in the committed colorectal table.
EXPECTED_DEFICIT = 0            # no taxon may select FEWER than were committed
EXPECTED_MEDIAN_RATIO = 1.11    # selected / committed, tolerance +- 0.06
EXPECTED_TAXA = 211

DRIFT = []


def load_exposure(path):
    """MiBioGen allHits file -> {taxon: [snp dicts]}."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bac = r["bac"].strip('"')
            out.setdefault(bac, []).append(dict(
                chrom=int(r["chr"]), bp=int(r["bp"]),
                rsid=r["rsID"].strip('"'),
                ref=r["ref.allele"].strip('"').upper(),
                eff=r["eff.allele"].strip('"').upper(),
                beta=float(r["beta"]), se=float(r["SE"]),
                p=float(r["P.weightedSumZ"]), n=float(r["N"])))
    return out


def select_instruments(hits, pthr=P_THRESHOLD, window=CLUMP_WINDOW,
                       drop_palindromic=DROP_PALINDROMIC):
    """p-threshold, palindromic removal, then greedy distance-based clumping."""
    keep = [h for h in hits if h["p"] < pthr]
    if drop_palindromic:
        keep = [h for h in keep
                if frozenset({h["ref"], h["eff"]}) not in PALINDROMIC]
    keep.sort(key=lambda h: h["p"])
    chosen = []
    for h in keep:
        if window and any(h["chrom"] == c["chrom"] and abs(h["bp"] - c["bp"]) < window
                          for c in chosen):
            continue
        chosen.append(h)
    return chosen


# ----------------------------- estimators ---------------------------------
def ivw(bx, by, sy):
    w = 1.0 / sy ** 2
    b = float(np.sum(w * bx * by) / np.sum(w * bx ** 2))
    se = float(math.sqrt(1.0 / np.sum(w * bx ** 2)))
    return b, se, 2 * stats.norm.sf(abs(b / se))


def mr_egger(bx, by, sy):
    w = 1.0 / sy ** 2
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    coef = cov @ (X.T @ W @ by)
    resid = by - X @ coef
    dof = max(len(bx) - 2, 1)
    sigma2 = float((resid ** 2 * w).sum() / dof)
    seb = np.sqrt(np.diag(cov) * sigma2)
    return (float(coef[1]), float(seb[1]),
            2 * stats.t.sf(abs(coef[1] / seb[1]), dof),
            float(coef[0]), 2 * stats.t.sf(abs(coef[0] / seb[0]), dof))


def weighted_median(bx, by, sy):
    ratio = by / bx
    w = (bx ** 2) / (sy ** 2)
    o = np.argsort(ratio)
    r, w = ratio[o], w[o]
    cw = np.cumsum(w) - 0.5 * w
    cw /= w.sum()
    k = np.searchsorted(cw, 0.5)
    k = min(max(k, 1), len(r) - 1)
    b = r[k - 1] + (r[k] - r[k - 1]) * (0.5 - cw[k - 1]) / (cw[k] - cw[k - 1])
    boot = []
    rng = np.random.default_rng(0)
    for _ in range(1000):
        bys = rng.normal(by, sy)
        rs = bys / bx
        o2 = np.argsort(rs)
        c2 = np.cumsum(w[o2]) - 0.5 * w[o2]
        c2 /= w[o2].sum()
        k2 = min(max(np.searchsorted(c2, 0.5), 1), len(rs) - 1)
        rr = rs[o2]
        boot.append(rr[k2 - 1] + (rr[k2] - rr[k2 - 1]) *
                    (0.5 - c2[k2 - 1]) / (c2[k2] - c2[k2 - 1]))
    se = float(np.std(boot))
    return float(b), se, 2 * stats.norm.sf(abs(b / se))


def cochran_q(bx, by, sy, b_ivw):
    w = 1.0 / sy ** 2
    q = float(np.sum(w * (by - b_ivw * bx) ** 2))
    return q, stats.chi2.sf(q, max(len(bx) - 1, 1))


def bh(p, q=0.05):
    p = np.asarray(p)
    o = np.argsort(p)
    m = p.size
    ok = p[o] <= q * np.arange(1, m + 1) / m
    out = np.zeros(m, bool)
    if ok.any():
        out[o[:np.max(np.where(ok)[0]) + 1]] = True
    return out


# --------------------------- stage 1 validation ---------------------------
def validate_instruments(exposure_path, committed_path):
    """Assert exposure selection is a superset of the committed instruments.

    The committed n_snp is the count AFTER harmonisation with the outcome
    GWAS, so exposure selection must yield at least as many for every taxon.
    A deficit anywhere means the selection parameters are wrong.
    """
    exposure = load_exposure(exposure_path)
    truth = {r["taxon"]: int(r["n_snp"])
             for r in csv.DictReader(open(committed_path), delimiter="\t")}

    if len(truth) != EXPECTED_TAXA:
        DRIFT.append(f"DRIFT taxa count: got {len(truth)}, expected {EXPECTED_TAXA}")
    missing = set(truth) - set(exposure)
    if missing:
        DRIFT.append(f"DRIFT {len(missing)} committed taxa absent from exposure file")

    deficits, ratios, exact = [], [], 0
    for tax, n in truth.items():
        if tax not in exposure:
            continue
        s = len(select_instruments(exposure[tax]))
        ratios.append(s / max(n, 1))
        if s < n:
            deficits.append((tax, s, n))
        elif s == n:
            exact += 1

    med = st.median(ratios)
    print(f"taxa checked                 : {len(ratios)}")
    print(f"selected >= committed        : {len(ratios) - len(deficits)}/{len(ratios)}")
    print(f"selected == committed        : {exact}/{len(ratios)}")
    print(f"median selected/committed    : {med:.2f}")
    print(f"parameters                   : p<{P_THRESHOLD:g}, "
          f"{CLUMP_WINDOW // 1000}kb window, palindromic "
          f"{'dropped' if DROP_PALINDROMIC else 'kept'}")

    if len(deficits) != EXPECTED_DEFICIT:
        DRIFT.append(f"DRIFT deficit count: got {len(deficits)}, "
                     f"expected {EXPECTED_DEFICIT}; e.g. {deficits[:3]}")
    if abs(med - EXPECTED_MEDIAN_RATIO) > 0.06:
        DRIFT.append(f"DRIFT median ratio: got {med:.3f}, "
                     f"expected {EXPECTED_MEDIAN_RATIO} +- 0.06")
    return len(ratios), exact, med


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    # Stage 2 is NOT here. See src/wp1_stage2_validate.py (CORRECTIONS C21).

    v = sub.add_parser("validate-instruments")
    v.add_argument("--exposure", required=True)
    v.add_argument("--committed", default="results/MR_results_colorectal.tsv")

    a = ap.parse_args()

    if a.cmd == "validate-instruments":
        for p in (a.exposure, a.committed):
            if not os.path.exists(p):
                print(f"missing input: {p}", file=sys.stderr)
                return 2
        validate_instruments(a.exposure, a.committed)
        if DRIFT:
            print("\n".join(DRIFT), file=sys.stderr)
            return 1
        print("\nStage 1 (exposure) validated against all committed instrument counts.")
        return 0

if __name__ == "__main__":
    sys.exit(main())
