#!/usr/bin/env python3
"""
RO-2026-008 | WP1 | Cross-cohort consistency and per-taxon triangulation
wp1_cross_consistency-1.0.0

Purpose
-------
Canonicalises an analysis that was previously run inline and never persisted
(reproducibility gap D3). Closes two further gaps:

  D2  The reported 61% cross-cohort directional agreement was tie-dependent.
      This module recovers the sign of every taxon, including those whose OR
      rounds to 1.00, so the agreement statistic no longer depends on a tie
      rule.

  D4  MR_results_colorectal.tsv stores the effect as a rounded display string
      "OR (LCI-UCI)" with no beta and no standard error. This module
      reconstructs b and se from the confidence interval and VALIDATES the
      reconstruction against the stored p-value before using it.

Method
------
For the FinnGen colorectal table, the log-scale effect is recovered as the
midpoint of the log confidence interval rather than as log(OR):

    b_hat  = (ln(LCI) + ln(UCI)) / 2
    se_hat = (ln(UCI) - ln(LCI)) / (2 * 1.959964)

The log-CI midpoint carries more information than the 2-decimal OR because it
is built from two independently rounded bounds. An OR displayed as 1.00 with
an asymmetric interval still yields a signed b_hat. Every reconstruction is
checked against the stored p_ivw; taxa failing tolerance are excluded and
reported, never silently used.

Outputs
-------
  cross_consistency.json         machine-readable results
  MR_cross_cohort_recon.tsv      reconstructed FinnGen effects with QC flags
  Figure7_cross_cohort.png       log-OR concordance scatter, 300 dpi

Usage
-----
    python3 wp1_cross_consistency.py [--outdir .] [--tol 0.15]

Exit codes
----------
    0  all reported values re-derived within tolerance
    1  a reported value failed to re-derive (DRIFT)
"""

import argparse
import json
import math
import sys

import numpy as np
from scipy import stats

Z = 1.959964  # two-sided 95% normal quantile

# Values asserted by this module. Any drift fails the run.
#
# n_common is 210, not 211: phylum.Cyanobacteria.id.1500 is excluded by
# reconstruction QC. It carries the smallest p in the FinnGen table
# (1.16e-4), and CI-based reconstruction loses precision in the tail, so its
# reconstructed p (6.8e-5) falls outside the log10 tolerance. It is excluded
# rather than used. This costs nothing analytically: that taxon was already
# diagnosed as a pleiotropic artifact (NOS2 instrument rs2314810) and its
# direction flipped on replication.
EXPECTED = {
    "n_common": 210,
    "n_qc_fail": 1,
    "agreement_recovered_pct": 61.4,   # tolerance +- 2.0 percentage points
    "alistipes_finngen_or": 0.95,      # displayed, 2dp; tolerance +- 0.015
    "alistipes_finngen_p": 0.5777,     # tolerance +- 0.01
    "alistipes_100k_or": 0.9334,       # tolerance +- 0.001
    "alistipes_100k_p": 0.0478,        # tolerance +- 0.001
    "alistipes_panc_or": 0.6823,       # tolerance +- 0.001
    "alistipes_panc_p": 0.0168,        # tolerance +- 0.001
}

DRIFT = []


def check(name, got, expected, tol):
    ok = abs(got - expected) <= tol
    if not ok:
        DRIFT.append(f"DRIFT {name}: got {got!r}, expected {expected!r} +- {tol}")
    return ok


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def parse_or_ci(text):
    """Parse 'OR (LCI-UCI)' into (or_disp, lci, uci)."""
    or_part, ci_part = text.split(" (")
    lci, uci = ci_part.rstrip(")").split("-")
    return float(or_part), float(lci), float(uci)


def reconstruct(lci, uci):
    """Recover log-scale effect and standard error from a 95% CI."""
    b = (math.log(lci) + math.log(uci)) / 2.0
    se = (math.log(uci) - math.log(lci)) / (2.0 * Z)
    return b, se


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="results",
                help="directory holding the MR_results_*.tsv tables")
    ap.add_argument("--figdir", default="figures",
                help="directory to write Figure7 into")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--tol", type=float, default=0.15,
                    help="max absolute difference in log10(p) for a "
                         "reconstruction to be accepted")
    args = ap.parse_args()

    # ---- Cohort 1: FinnGen colorectal (reconstruction required) -----------
    fin = {}
    recon_rows = []
    n_fail = 0
    for r in read_tsv(f"{args.datadir}/MR_results_colorectal.tsv"):
        taxon = r["taxon"]
        or_disp, lci, uci = parse_or_ci(r["OR (95% CI)"])
        p_stored = float(r["p_ivw"])
        b, se = reconstruct(lci, uci)
        p_recon = 2.0 * stats.norm.sf(abs(b) / se) if se > 0 else float("nan")

        # QC: reconstruction must reproduce the stored p-value on log10 scale
        d = abs(math.log10(max(p_recon, 1e-300)) - math.log10(max(p_stored, 1e-300)))
        ok = d <= args.tol
        if not ok:
            n_fail += 1
        else:
            fin[taxon] = {"b": b, "se": se, "or": math.exp(b),
                          "p_stored": p_stored, "n_snp": int(r["n_snp"])}
        recon_rows.append({
            "taxon": taxon, "n_snp": r["n_snp"], "or_displayed": or_disp,
            "lci": lci, "uci": uci, "b_recon": round(b, 6),
            "se_recon": round(se, 6), "p_stored": p_stored,
            "p_recon": p_recon, "log10p_delta": round(d, 4),
            "qc": "PASS" if ok else "FAIL",
        })

    # ---- Cohort 2: Fernandez-Rozadilla 100k colorectal (b, se stored) -----
    big = {}
    for r in read_tsv(f"{args.datadir}/MR_results_colorectal_100k.tsv"):
        big[r["taxon"]] = {"b": float(r["b"]), "se": float(r["se"]),
                           "or": float(r["OR"]), "p": float(r["p"]),
                           "p_wm": float(r["p_wm"]), "n_snp": int(r["n"])}

    # ---- Cohort 3: FinnGen pancreatic ------------------------------------
    panc = {}
    for r in read_tsv(f"{args.datadir}/MR_results_pancreatic.tsv"):
        panc[r["taxon"]] = {"b": float(r["b_ivw"]), "se": float(r["se_ivw"]),
                            "or": float(r["or_ivw"]), "p": float(r["p_ivw"]),
                            "p_wm": float(r["p_wmedian"]), "n_snp": int(r["n_snp"])}

    # ---- Directional agreement, colorectal cohorts ------------------------
    common = sorted(set(fin) & set(big))
    n = len(common)
    agree = [t for t in common if fin[t]["b"] * big[t]["b"] > 0]
    k = len(agree)
    pct = 100.0 * k / n
    p_binom = stats.binomtest(k, n, 0.5).pvalue

    # Sensitivity: the two naive tie rules previously in conflict
    naive_or = {t: parse_or_ci(r["OR (95% CI)"])[0]
                for r in read_tsv(f"{args.datadir}/MR_results_colorectal.tsv")
                for t in [r["taxon"]]}
    ties = [t for t in common if naive_or[t] == 1.00 or big[t]["or"] == 1.00]
    k_ties_as_disagree = sum(1 for t in common
                             if (naive_or[t] - 1) * (big[t]["or"] - 1) > 0)
    nontie = [t for t in common if t not in ties]
    k_ties_dropped = sum(1 for t in nontie
                         if (naive_or[t] - 1) * (big[t]["or"] - 1) > 0)

    # ---- Per-taxon triangulation across all three tests -------------------
    tri_common = sorted(set(fin) & set(big) & set(panc))
    triangulated = []
    for t in tri_common:
        bs = [fin[t]["b"], big[t]["b"], panc[t]["b"]]
        ps = [fin[t]["p_stored"], big[t]["p"], panc[t]["p"]]
        same_dir = all(x < 0 for x in bs) or all(x > 0 for x in bs)
        if same_dir:
            triangulated.append({
                "taxon": t,
                "direction": "protective" if bs[0] < 0 else "risk",
                "or_crc_finngen": round(fin[t]["or"], 3),
                "p_crc_finngen": fin[t]["p_stored"],
                "or_crc_100k": round(big[t]["or"], 3),
                "p_crc_100k": big[t]["p"],
                "p_wm_crc_100k": big[t]["p_wm"],
                "or_panc": round(panc[t]["or"], 3),
                "p_panc": panc[t]["p"],
                "p_wm_panc": panc[t]["p_wm"],
                "n_nominal": sum(1 for p in ps if p < 0.05),
            })
    triangulated.sort(key=lambda d: (-d["n_nominal"], min(
        d["p_crc_finngen"], d["p_crc_100k"], d["p_panc"])))

    # ---- Assertions -------------------------------------------------------
    check("n_common", n, EXPECTED["n_common"], 0)
    check("n_qc_fail", n_fail, EXPECTED["n_qc_fail"], 0)
    check("agreement_recovered_pct", round(pct, 1),
          EXPECTED["agreement_recovered_pct"], 2.0)
    ali = "genus.Alistipes.id.968"
    # Compare the unrounded reconstruction against the 2dp displayed value.
    # Rounding before comparison would fail spuriously on the half-ulp.
    check("alistipes_finngen_or", fin[ali]["or"],
          EXPECTED["alistipes_finngen_or"], 0.015)
    check("alistipes_finngen_p", fin[ali]["p_stored"],
          EXPECTED["alistipes_finngen_p"], 0.01)
    check("alistipes_100k_or", big[ali]["or"], EXPECTED["alistipes_100k_or"], 0.001)
    check("alistipes_100k_p", big[ali]["p"], EXPECTED["alistipes_100k_p"], 0.001)
    check("alistipes_panc_or", panc[ali]["or"], EXPECTED["alistipes_panc_or"], 0.001)
    check("alistipes_panc_p", panc[ali]["p"], EXPECTED["alistipes_panc_p"], 0.001)

    # ---- Emit -------------------------------------------------------------
    out = {
        "module": "wp1_cross_consistency-1.0.0",
        "reconstruction": {
            "method": "log-CI midpoint; se from CI width / (2 * 1.959964)",
            "n_taxa": len(recon_rows),
            "n_qc_pass": len(fin),
            "n_qc_fail": n_fail,
            "tolerance_log10p": args.tol,
        },
        "directional_agreement": {
            "n_common": n,
            "n_agree": k,
            "pct": round(pct, 1),
            "binomial_p": p_binom,
            "tie_sensitivity": {
                "n_ties_under_naive_rule": len(ties),
                "naive_ties_counted_as_disagreement": {
                    "pct": round(100.0 * k_ties_as_disagree / n, 1),
                    "p": stats.binomtest(k_ties_as_disagree, n, 0.5).pvalue,
                },
                "naive_ties_dropped": {
                    "pct": round(100.0 * k_ties_dropped / len(nontie), 1),
                    "p": stats.binomtest(k_ties_dropped, len(nontie), 0.5).pvalue,
                },
            },
        },
        "triangulation": {
            "n_taxa_all_three": len(tri_common),
            "n_consistent_direction": len(triangulated),
            "expected_under_independence": round(len(tri_common) * 0.25, 1),
            "top": triangulated[:10],
        },
        "drift": DRIFT,
    }
    with open(f"{args.outdir}/cross_consistency.json", "w") as fh:
        json.dump(out, fh, indent=2)

    cols = ["taxon", "n_snp", "or_displayed", "lci", "uci", "b_recon",
            "se_recon", "p_stored", "p_recon", "log10p_delta", "qc"]
    with open(f"{args.outdir}/MR_cross_cohort_recon.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recon_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    pcts = {
        "recovered": pct,
        "ties_dropped": 100.0 * k_ties_dropped / len(nontie),
        "ties_disagree": 100.0 * k_ties_as_disagree / n,
    }
    make_figure(common, fin, big, agree, pcts, args.figdir)

    # ---- Report -----------------------------------------------------------
    print(f"Reconstruction QC: {len(fin)}/{len(recon_rows)} pass, {n_fail} fail")
    print(f"Directional agreement (sign-recovered): {k}/{n} = {pct:.1f}%, "
          f"binomial p = {p_binom:.4g}")
    print(f"  naive rule, ties as disagreement: "
          f"{out['directional_agreement']['tie_sensitivity']['naive_ties_counted_as_disagreement']['pct']}%")
    print(f"  naive rule, ties dropped:         "
          f"{out['directional_agreement']['tie_sensitivity']['naive_ties_dropped']['pct']}%")
    print(f"Consistent-direction taxa across all three tests: "
          f"{len(triangulated)}/{len(tri_common)} "
          f"(expected under independence: {len(tri_common) * 0.25:.1f})")
    print("\nTop triangulated taxa:")
    for d in triangulated[:5]:
        print(f"  {d['taxon']:<45} {d['direction']:<10} "
              f"OR {d['or_crc_finngen']}/{d['or_crc_100k']}/{d['or_panc']}  "
              f"nominal in {d['n_nominal']}/3")

    if DRIFT:
        print("\n".join(DRIFT), file=sys.stderr)
        return 1
    print("\nAll asserted values re-derived within tolerance.")
    return 0


def make_figure(common, fin, big, agree, pcts, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.array([fin[t]["b"] for t in common])
    y = np.array([big[t]["b"] for t in common])
    ok = np.array([t in set(agree) for t in common])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel a: concordance scatter. Okabe-Ito, colorblind safe.
    ax[0].axhline(0, color="0.75", lw=0.8, zorder=0)
    ax[0].axvline(0, color="0.75", lw=0.8, zorder=0)
    ax[0].scatter(x[ok], y[ok], s=22, c="#0072B2", alpha=0.75,
                  label=f"same direction (n={ok.sum()})", edgecolors="none")
    ax[0].scatter(x[~ok], y[~ok], s=22, c="#D55E00", alpha=0.75, marker="^",
                  label=f"opposite direction (n={(~ok).sum()})", edgecolors="none")
    ali = "genus.Alistipes.id.968"
    if ali in fin and ali in big:
        ax[0].scatter([fin[ali]["b"]], [big[ali]["b"]], s=90,
                      facecolors="none", edgecolors="#000000", lw=1.4, zorder=5)
        ax[0].annotate("Alistipes", (fin[ali]["b"], big[ali]["b"]),
                       textcoords="offset points", xytext=(8, -12), fontsize=9)
    ax[0].set_xlabel("log OR, FinnGen R12 (11,790 cases)")
    ax[0].set_ylabel("log OR, Fernandez-Rozadilla\n(100,204 cases)")
    ax[0].set_title(f"a  Cross-cohort concordance, {len(common)} gut taxa",
                    loc="left", fontsize=11)
    ax[0].legend(frameon=False, fontsize=8, loc="upper left")

    # Panel b: agreement against the null under all three tie rules.
    # The point is robustness: the result does not depend on the rule.
    labels = ["sign-\nrecovered", "naive,\nties dropped", "naive,\nties as disagree"]
    vals = [pcts["recovered"], pcts["ties_dropped"], pcts["ties_disagree"]]
    colors = ["#0072B2", "#56B4E9", "#56B4E9"]
    ax[1].bar(range(3), vals, width=0.62, color=colors)
    ax[1].axhline(50, color="#D55E00", ls="--", lw=1.4)
    ax[1].text(2.52, 50.5, "expected\nunder noise", color="#D55E00", fontsize=8,
               ha="right", va="bottom")
    ax[1].set_ylim(45, 68)
    ax[1].set_xticks(range(3))
    ax[1].set_xticklabels(labels, fontsize=8.5)
    ax[1].set_ylabel("directional agreement (%)")
    ax[1].set_title("b  Excess over noise is robust to the tie rule",
                    loc="left", fontsize=11)
    for i, v in enumerate(vals):
        ax[1].text(i, v + 0.5, f"{v:.1f}%", ha="center", fontsize=9)

    for a in ax:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(f"{outdir}/Figure7_cross_cohort.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
