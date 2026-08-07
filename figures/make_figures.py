#!/usr/bin/env python3
"""
RO-2026-008 | Figure regeneration from committed results
make_figures-1.0.0

Closes limitation L1 for Figures 1, 2, 3 and 6. Those figures were originally
produced interactively and their plotting code was lost; only the PNGs
survived. The underlying results were committed, so the figures are rebuilt
here from results/ at any resolution.

Figures 4 and 5 still require the TCMA download described in docs/RUNBOOK.md
and are not covered by this module.

Faithfulness note: this is a RECONSTRUCTION of the plotting code, not a
recovery of it. The data plotted is identical and assertion-checked against
docs/CANONICAL_RESULTS.md. Cosmetic details may differ from the originals.

Usage
-----
    python3 figures/make_figures.py [--dpi 600] [--outdir figures]

Exit codes
----------
    0  all asserted values re-derived within tolerance
    1  a value failed to re-derive (DRIFT)
"""

import argparse
import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito, colorblind safe
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#999999"
CHANCE = 1.0 / 6.0          # six cancer types
FLOOR = 0.8                 # detection floor, log units

DRIFT = []


def check(name, got, expected, tol):
    if abs(got - expected) > tol:
        DRIFT.append(f"DRIFT {name}: got {got!r}, expected {expected!r} +- {tol}")


def load_audit(d):
    rows = []
    for f in ("audit_seeds_0_4.json", "audit_seeds_5_9.json"):
        rows += json.load(open(os.path.join(d, f)))["rows"]
    return rows


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# --------------------------------------------------------------------------
def figure1(rows, out, dpi):
    """Which tests catch fabricated signal, and which metrics separate."""
    tests = [("t1", "T1\nno-info rate"), ("t2", "T2\npermutation"),
             ("t3", "T3\nconfounder"), ("t5a", "T5a\nwithin-batch")]
    A = [r for r in rows if r["scenario"] == "A"]
    B = [r for r in rows if r["scenario"] == "B"]
    passA = [sum(bool(r[t]) for r in A) for t, _ in tests]
    passB = [sum(bool(r[t]) for r in B) for t, _ in tests]

    check("F1 t1 pass A", passA[0], 10, 0); check("F1 t1 pass B", passB[0], 10, 0)
    check("F1 t2 pass A", passA[1], 10, 0); check("F1 t2 pass B", passB[1], 10, 0)
    check("F1 t3 pass A", passA[2], 10, 0); check("F1 t3 pass B", passB[2], 0, 0)
    check("F1 t5a pass A", passA[3], 10, 0); check("F1 t5a pass B", passB[3], 0, 0)

    fig, ax = plt.subplots(1, 2, figsize=(9.75, 5.0))
    fig.suptitle("Figure 1  Conventional significance tests pass on data with no biology",
                 fontsize=11, x=0.01, ha="left", y=0.985)

    x = np.arange(len(tests)); w = 0.38
    ax[0].bar(x - w / 2, passA, w, color=BLUE, label="real signal (cohort A)")
    ax[0].bar(x + w / 2, passB, w, color=ORANGE, label="zero signal (cohort B)")
    for xi, (a, b) in enumerate(zip(passA, passB)):
        ax[0].text(xi - w / 2, a + 0.18, str(a), ha="center", fontsize=8.5)
        ax[0].text(xi + w / 2, b + 0.18, str(b), ha="center", fontsize=8.5)
    ax[0].annotate("pass on\nNOISE", xy=(1 + w / 2, 10), xytext=(1.60, 8.7),
                   fontsize=8.5, color=ORANGE, ha="left",
                   arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
    ax[0].set_xticks(x); ax[0].set_xticklabels([l for _, l in tests], fontsize=8.5)
    ax[0].set_ylabel("pass rate (of 10 seeds)")
    ax[0].set_ylim(0, 11.4)
    ax[0].set_title("a  Which tests catch fabricated signal", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.34))
    style(ax[0])

    data = [[r["t3_margin"] for r in A], [r["t3_margin"] for r in B],
            [r["t5a_delta"] for r in A], [r["t5a_delta"] for r in B]]
    pos = [0.8, 1.3, 2.4, 2.9]
    bp = ax[1].boxplot(data, positions=pos, widths=0.36, patch_artist=True,
                       medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], [BLUE, ORANGE, BLUE, ORANGE]):
        patch.set_facecolor(c); patch.set_alpha(0.80); patch.set_edgecolor("none")
    rng = np.random.default_rng(0)
    for series, p in zip(data, pos):
        ax[1].plot(p + rng.uniform(-0.09, 0.09, len(series)), series, ".",
                   color="#333333", ms=3.6, alpha=0.75, zorder=3)
    ax[1].axhline(0, color=GREY, ls="--", lw=1.1)
    ax[1].set_xticks([1.05, 2.65]); ax[1].set_xticklabels(["T3 margin", "T5a delta"])
    ax[1].set_ylabel("metric value")
    ax[1].set_title("b  Discriminating metrics separate cleanly", loc="left", fontsize=10)
    ax[1].text(0.8, max(data[0]) + 0.02, "real", ha="center", fontsize=8, color=BLUE)
    ax[1].text(2.9, min(data[3]) - 0.05, "zero signal", ha="center", fontsize=8, color=ORANGE)
    style(ax[1])

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(out, "Figure1_audit_validation.png"), dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------
def figure2(sweep, out, dpi):
    """Confounding manufactures accuracy; the within-batch test resists it."""
    conf = sorted({r["confound"] for r in sweep})
    def series(sc, key):
        return [next(r[key] for r in sweep if r["scenario"] == sc and r["confound"] == c)
                for c in conf]

    accB = series("B", "acc_mean")
    check("F2 zero-signal acc at min confounding", accB[0], 0.169, 0.002)
    check("F2 zero-signal acc at max confounding", accB[-1], 0.409, 0.002)
    ratio = accB[-1] / CHANCE

    fig, ax = plt.subplots(1, 2, figsize=(9.75, 5.0))
    fig.suptitle("Figure 2  Batch confounding alone produces false signal; the within-batch test resists it",
                 fontsize=11, x=0.01, ha="left", y=0.985)

    ax[0].errorbar(conf, series("A", "acc_mean"), yerr=series("A", "acc_sd"),
                   marker="o", color=BLUE, capsize=3, label="real signal")
    ax[0].errorbar(conf, accB, yerr=series("B", "acc_sd"),
                   marker="o", color=ORANGE, capsize=3, label="zero signal")
    ax[0].axhline(CHANCE, color=GREY, ls="--", lw=1.2)
    ax[0].text(conf[-1], CHANCE + 0.014, f"chance ({CHANCE:.3f})", fontsize=8,
               color=GREY, ha="right")
    ax[0].annotate(f"{ratio:.1f}x chance\nfrom ZERO biology",
                   xy=(conf[-1], accB[-1]), xytext=(conf[-1] - 0.30, accB[-1] + 0.06),
                   fontsize=8.5, color=ORANGE,
                   arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
    ax[0].set_xlabel("batch-outcome confounding strength")
    ax[0].set_ylabel("classification accuracy")
    ax[0].set_title("a  Confounding manufactures accuracy from noise", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    style(ax[0])

    ax[1].axhspan(-0.05, 0.05, color=GREY, alpha=0.18)
    ax[1].errorbar(conf, series("A", "t5a_delta_mean"), yerr=series("A", "t5a_delta_sd"),
                   marker="o", color=BLUE, capsize=3, label="real signal")
    ax[1].errorbar(conf, series("B", "t5a_delta_mean"), yerr=series("B", "t5a_delta_sd"),
                   marker="o", color=ORANGE, capsize=3, label="zero signal")
    ax[1].axhline(0, color=GREY, ls="--", lw=1.1)
    ax[1].set_xlabel("batch-outcome confounding strength")
    ax[1].set_ylabel("within-batch delta (T5a)")
    ax[1].set_title("b  Within-batch test is immune to confounding", loc="left", fontsize=10)
    ax[1].legend(frameon=False, fontsize=8, loc="upper left")
    style(ax[1])

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(out, "Figure2_confounding_sweep.png"), dpi=dpi)
    plt.close(fig)
    return ratio


# --------------------------------------------------------------------------
def figure3(floor, out, dpi):
    """Detection floor by true effect size."""
    e = [r["effect"] for r in floor]
    sens = [r["sensitivity_mean"] for r in floor]
    sd = [r["sensitivity_sd"] for r in floor]
    fd = [r["false_disc_mean"] for r in floor]

    check("F3 sensitivity at 0.4", sens[e.index(0.4)], 0.235, 0.005)
    check("F3 sensitivity at 0.8", sens[e.index(0.8)], 0.98, 0.02)

    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    fig.suptitle(f"Figure 3  Detection floor \u2248 {FLOOR} log units",
                 fontsize=11, x=0.02, ha="left", y=0.985)
    ax.errorbar(e, sens, yerr=sd, marker="o", color=BLUE, capsize=3, label="sensitivity")
    ax.plot(e, fd, marker="s", ls="--", color=ORANGE, label="false discovery")
    ax.axhline(0.05, color=GREY, ls="--", lw=1.0)
    ax.text(0.105, 0.062, "nominal FDR 0.05", fontsize=7.5, color=GREY)
    ax.axvline(FLOOR, color=GREEN, ls="--", lw=1.2)
    ax.axhline(0.8, color=GREEN, ls=":", lw=1.1)
    ax.text(FLOOR + 0.02, 0.35, f"{FLOOR} sensitivity\nthreshold", fontsize=8, color=GREEN)
    ax.annotate("below floor:\nreal signal missed", xy=(0.235, 0.14),
                xytext=(0.30, 0.55), fontsize=8, color="#333333",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0))
    ax.set_xlabel("true effect size (log units)")
    ax.set_ylabel("rate")
    ax.set_ylim(-0.04, 1.08)
    ax.legend(frameon=False, fontsize=8.5, loc="center right")
    style(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(out, "Figure3_detection_floor.png"), dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------
def figure6(path, out, dpi, n_show=12):
    """Forest plot of the FinnGen colorectal MR screen."""
    rank = {"phylum": "p", "class": "c", "order": "o",
            "family": "f", "genus": "g", "species": "s"}
    recs = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            orv, ci = r["OR (95% CI)"].split(" (")
            lci, uci = [float(x) for x in ci.rstrip(")").split("-")]
            parts = r["taxon"].split(".")
            lvl = rank.get(parts[0], "?")
            name = parts[1] if len(parts) > 1 else r["taxon"]
            fdr = float(r["fdr_ivw"])
            pleio = float(r["egger_intercept_p"]) < 0.05
            het = float(r["q_p"]) < 0.05
            flag = pleio or het
            tag = ("FDR" if fdr < 0.05 else
                   ",".join([t for t, on in (("pleio", pleio), ("het", het)) if on])
                   or "clean")
            recs.append(dict(label=f"{lvl}_{name} (n={r['n_snp']}, {tag})",
                             orv=float(orv), lci=lci, uci=uci,
                             p=float(r["p_ivw"]), fdr=fdr, flag=flag))
    recs.sort(key=lambda d: d["p"])
    recs = recs[:n_show][::-1]

    fig, ax = plt.subplots(figsize=(9.55, 5.0))
    fig.suptitle("Figure 6  Two-sample MR: gut taxa \u2192 colorectal cancer (FinnGen, 11,790 cases)",
                 fontsize=11, x=0.01, ha="left", y=0.985)
    for i, d in enumerate(recs):
        c = ORANGE if d["fdr"] < 0.05 else (GREY if d["flag"] else BLUE)
        ax.plot([d["lci"], d["uci"]], [i, i], color=c, lw=1.8, solid_capstyle="round")
        ax.plot([d["orv"]], [i], "o", color=c, ms=7)
    ax.axvline(1.0, color="black", ls="--", lw=1.1)
    ax.set_yticks(range(len(recs)))
    ax.set_yticklabels([d["label"] for d in recs], fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks([0.6, 0.8, 1.0, 1.25, 1.5, 2.0])
    ax.set_xticklabels(["0.6", "0.8", "1.0", "1.25", "1.5", "2.0"])
    ax.set_xlabel("OR for colorectal cancer per SD increase in taxon abundance (IVW, 95% CI)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=ORANGE, lw=2, label="survives FDR"),
                       Line2D([], [], color=BLUE, lw=2, label="clean, suggestive only"),
                       Line2D([], [], color=GREY, lw=2, label="pleiotropy/heterogeneity flag")],
              frameon=False, fontsize=7.5, loc="lower right")
    style(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(out, "Figure6_MR_colorectal.png"), dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--datadir", default="results")
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = load_audit(a.datadir)
    sweep = json.load(open(os.path.join(a.datadir, "sweep_3seeds.json")))["rows"]
    floor = json.load(open(os.path.join(a.datadir, "floor_10seeds.json")))["rows"]

    figure1(rows, a.outdir, a.dpi)
    ratio = figure2(sweep, a.outdir, a.dpi)
    figure3(floor, a.outdir, a.dpi)
    figure6(os.path.join(a.datadir, "MR_results_colorectal.tsv"), a.outdir, a.dpi)

    print(f"Regenerated Figures 1, 2, 3, 6 at {a.dpi} dpi into {a.outdir}/")
    print(f"Confounding ratio at max sweep (0.95): {ratio:.2f}x chance")
    if DRIFT:
        print("\n".join(DRIFT), file=sys.stderr)
        return 1
    print("All asserted values re-derived within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
