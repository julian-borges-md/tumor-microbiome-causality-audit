#!/usr/bin/env python3
"""
RO-2026-008 | WP0 | Taxonomic redundancy sweep
wp0_nesting_sweep-1.0.0

Canonicalises nest_sweep.py, which was exploratory: it pulled its dependencies
in with `exec(open('wp0_nesting.py').read().split('if __name__')[0])`, emitted
only formatted text, and asserted nothing. MANIFEST.md flagged it as "not yet
canonicalised" and it was the last thing blocking a generator for Figure 5.

This module imports its dependencies properly, emits machine-readable results,
and fails on drift.

Measures two things on real TCMA data:

  A. Discovery inflation. How many "findings" a single analysis reports as the
     rule for collapsing nested taxonomic ranks is varied, from none (report
     every rank separately) to the pure ancestor rule (collapse any lineage
     that is a prefix of another).

  B. Structural redundancy of the feature space itself, independent of any
     result.

Requires the TCMA download described in docs/RUNBOOK.md.

    export TCMA_DIR=/path/to/tcma
    python3 src/wp0_nesting_sweep.py [--out results/nesting_sweep.json]

Exit codes
----------
    0  all asserted values re-derived
    1  a value failed to re-derive (DRIFT)
    2  TCMA_DIR not set
"""

import argparse
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wp0_nesting as nst                       # noqa: E402

R_THRESHOLDS = [0.99, 0.95, 0.90, 0.80]

# Values asserted by this module, from docs/findings/WP0_NESTING_FINDINGS.md
EXPECTED_DISCOVERY = {
    "STAD paired":   {"reported": 8,  "r>0.99": 3,  "r>0.95": 3,  "r>0.90": 2,  "r>0.80": 2,  "ancestor": 2},
    "HNSC unpaired": {"reported": 51, "r>0.99": 44, "r>0.95": 43, "r>0.90": 42, "r>0.80": 31, "ancestor": 9},
}
EXPECTED_FEATURES = {
    "COAD": {"nominal": 500, "r>0.99": 418, "r>0.95": 388, "r>0.90": 368},
    "STAD": {"nominal": 452, "r>0.99": 359, "r>0.95": 329, "r>0.90": 306},
    "HNSC": {"nominal": 500, "r>0.99": 402, "r>0.95": 384, "r>0.90": 366},
    "ESCA": {"nominal": 431, "r>0.99": 310, "r>0.95": 244, "r>0.90": 202},
}

DRIFT = []


def check(name, got, expected):
    if got != expected:
        DRIFT.append(f"DRIFT {name}: got {got}, expected {expected}")


def sig_set(X, meta, proj, paired=True):
    """Significant taxa at BH-FDR<0.05, plus the matrix they were tested on."""
    if paired:
        T, N, cols = nst.paired_matrix(X, meta, proj)
        if T is None or T.shape[0] < 5:
            return None, None, None
        nz = (T.sum(0) + N.sum(0)) > 0
        Ts, Ns, ids = T[:, nz], N[:, nz], cols[nz]
        pv = np.ones(Ts.shape[1])
        for j in range(Ts.shape[1]):
            d = Ts[:, j] - Ns[:, j]
            if np.count_nonzero(d) < 3:
                continue
            try:
                pv[j] = stats.wilcoxon(d, zero_method="wilcox").pvalue
            except Exception:
                pv[j] = 1.0
        M = np.vstack([Ts, Ns])
    else:
        m = meta["project"] == proj
        sub, subm = X.loc[m], meta.loc[m]
        tum = subm["sample_type"].str.contains("Tumor", na=False)
        nor = subm["sample_type"].str.contains("Normal", na=False)
        sel = tum | nor
        Xv = np.log1p(sub.loc[sel].values.astype(float))
        kc = Xv.sum(0) > 0
        Xv, ids = Xv[:, kc], np.array(sub.columns)[kc]
        y = tum[sel].astype(int).values
        pv = np.array([stats.mannwhitneyu(Xv[y == 1, j], Xv[y == 0, j],
                                          alternative="two-sided").pvalue
                       for j in range(Xv.shape[1])])
        M = Xv
    pv = np.nan_to_num(pv, nan=1.0)
    hit = nst.bh(pv, 0.05)
    return list(ids[hit]), M, list(ids)


def ancestor_clades(hids, tax):
    """Collapse any taxon whose lineage string is a prefix of another's."""
    lin = {i: (tax.loc[float(i), "lineage"] if float(i) in tax.index else str(i))
           for i in hids}
    clades = []
    for i in sorted(hids, key=lambda z: len(lin[z])):
        placed = False
        for c in clades:
            if lin[i].startswith(lin[c[0]]) or lin[c[0]].startswith(lin[i]):
                c.append(i)
                placed = True
                break
        if not placed:
            clades.append([i])
    return clades


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/nesting_sweep.json")
    a = ap.parse_args()
    if not os.environ.get("TCMA_DIR"):
        print("TCMA_DIR not set. See docs/RUNBOOK.md.", file=sys.stderr)
        return 2

    X, meta = nst.load()
    tax = nst.tax_table()
    out = {"module": "wp0_nesting_sweep-1.0.0", "discovery": {}, "features": {}}

    # ---- A. discovery inflation ------------------------------------------
    print("=== A. DISCOVERY COUNT vs REDUNDANCY RULE ===")
    hdr = f"{'analysis':<16}{'reported':>9}" + \
          "".join(f"{f'r>{t}':>8}" for t in R_THRESHOLDS) + f"{'ancestor':>10}{'fold':>7}"
    print(hdr)
    for label, proj, paired in [("STAD paired", "STAD", True),
                                ("HNSC unpaired", "HNSC", False)]:
        hids, M, allids = sig_set(X, meta, proj, paired)
        if not hids:
            continue
        Mh = M[:, [allids.index(h) for h in hids]]
        rec = {"reported": len(hids)}
        for t in R_THRESHOLDS:
            g, _ = nst.redundancy_graph(hids, tax, Mh, rthresh=t)
            rec[f"r>{t:.2f}"] = len(g)
        rec["ancestor"] = len(ancestor_clades(hids, tax))
        rec["fold"] = round(rec["reported"] / max(rec["ancestor"], 1), 2)
        out["discovery"][label] = rec
        for k, v in EXPECTED_DISCOVERY[label].items():
            check(f"{label} {k}", rec[k], v)
        print(f"{label:<16}{rec['reported']:>9}" +
              "".join(f"{rec[f'r>{t:.2f}']:>8}" for t in R_THRESHOLDS) +
              f"{rec['ancestor']:>10}{rec['fold']:>6.1f}x")

    # ---- B. feature-space redundancy -------------------------------------
    print("\n=== B. FEATURE-SPACE REDUNDANCY (paired cohorts, first 500 taxa) ===")
    print(f"{'cancer':>7}{'nominal':>9}" + "".join(f"{f'r>{t}':>8}" for t in R_THRESHOLDS[:3]))
    for proj in ["COAD", "STAD", "HNSC", "ESCA"]:
        T, N, cols = nst.paired_matrix(X, meta, proj)
        if T is None:
            continue
        nz = (T.sum(0) + N.sum(0)) > 0
        ids = list(cols[nz])
        M = np.vstack([T[:, nz], N[:, nz]])
        if len(ids) > 500:
            ids, M = ids[:500], M[:, :500]
        rec = {"nominal": len(ids)}
        for t in R_THRESHOLDS[:3]:
            g, _ = nst.redundancy_graph(ids, tax, M, rthresh=t)
            rec[f"r>{t:.2f}"] = len(g)
        rec["redundant_pct_at_090"] = round(100 * (1 - rec["r>0.90"] / rec["nominal"]))
        out["features"][proj] = rec
        for k, v in EXPECTED_FEATURES[proj].items():
            check(f"{proj} {k}", rec[k], v)
        print(f"{proj:>7}{rec['nominal']:>9}" +
              "".join(f"{rec[f'r>{t:.2f}']:>8}" for t in R_THRESHOLDS[:3]))

    out["drift"] = DRIFT
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {a.out}")

    if DRIFT:
        print("\n".join(DRIFT), file=sys.stderr)
        return 1
    print("All asserted values re-derived.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
