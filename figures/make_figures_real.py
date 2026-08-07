#!/usr/bin/env python3
"""
RO-2026-008 | Real-data figure regeneration
make_figures_real-1.0.0

Rebuilds Figures 4 and 5 from The Cancer Microbiome Atlas, closing the last
of limitation L1. Requires the TCMA download described in
docs/RUNBOOK.md.

    export TCMA_DIR=/path/to/tcma
    python3 figures/make_figures_real.py --dpi 600 --outdir figures

Faithfulness note. This reconstructs plotting code that was lost; the analysis
code it calls is the committed pipeline, unchanged. Panel b (the paired
H. pylori series) re-derives exactly against the canonical value. Panel a
plots cross-validated ROC AUC, which was never entered into
docs/CANONICAL_RESULTS.md and moves by up to 0.03 across scikit-learn
versions. The claim it supports, weak within-tissue discrimination in all four
evaluable cancers, is unaffected: every AUC sits between 0.50 and 0.70 under
both versions. See docs/CORRECTIONS.md C9.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import wp0_tumor_vs_normal as tvn          # noqa: E402
import wp0_paired as pw                    # noqa: E402
import wp0_nesting_sweep as nsw            # noqa: E402

BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#999999"
PROJECTS = ["COAD", "ESCA", "STAD", "HNSC"]
DRIFT = []


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_a(ax, X, meta):
    aucs, ns = [], []
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    for proj in PROJECTS:
        m = meta["project"] == proj
        sub, subm = X.loc[m], meta.loc[m]
        st = subm["sample_type"]
        tum = st.str.contains("Tumor", na=False)
        nor = st.str.contains("Normal", na=False)
        sel = tum | nor
        Xv = np.log1p(sub.loc[sel].values.astype(float))
        Xv = Xv[:, Xv.sum(0) > 0]
        y = tum[sel].astype(int).values
        clf = RandomForestClassifier(n_estimators=tvn.N_TREES, n_jobs=-1,
                                     random_state=0, min_samples_leaf=2)
        aucs.append(cross_val_score(clf, Xv, y, cv=cv, scoring="roc_auc", n_jobs=1).mean())
        ns.append((int(y.sum()), int((1 - y).sum())))

    for p, (t, n), exp in zip(PROJECTS, ns, [(125, 21), (62, 22), (128, 39), (157, 22)]):
        if (t, n) != exp:
            DRIFT.append(f"DRIFT {p} sample counts: got {(t, n)}, expected {exp}")

    colors = [GREY if a < 0.60 else GREEN for a in aucs]
    ax.bar(range(len(PROJECTS)), aucs, width=0.62, color=colors)
    for i, (a, (t, n)) in enumerate(zip(aucs, ns)):
        ax.text(i, a + 0.023, f"{a:.2f}", ha="center", fontsize=9)
        ax.text(i, a + 0.006, f"n={t}v{n}", ha="center", fontsize=7.5, color="#444444")
    ax.axhline(0.5, color="black", ls="--", lw=1.2)
    ax.text(len(PROJECTS) - 0.55, 0.507, "no signal (0.5)", fontsize=8, ha="right")
    ax.set_xticks(range(len(PROJECTS)))
    ax.set_xticklabels(PROJECTS)
    ax.set_ylabel("tumor vs adjacent-normal AUC")
    ax.set_ylim(0.45, max(aucs) + 0.075)
    ax.set_title("a  Within-tissue tumor signal is weak", loc="left", fontsize=10)
    style(ax)
    return aucs, ns


def panel_b(ax, X, meta):
    tn = pw.taxnames()
    name2id = {v: k for k, v in tn.items() if isinstance(v, str)}
    T, N, cols = pw.build_pairs(X, meta, "STAD")
    j = [float(c) for c in cols].index(name2id["Helicobacter pylori"])
    tum, nor = T[:, j], N[:, j]
    d = tum - nor
    diff = d.mean()
    p = stats.wilcoxon(tum, nor).pvalue

    if abs(diff - (-0.994)) > 0.002:
        DRIFT.append(f"DRIFT H. pylori paired diff: got {diff:.4f}, expected -0.994")
    if T.shape[0] != 39:
        DRIFT.append(f"DRIFT STAD pairs: got {T.shape[0]}, expected 39")

    for a, b in zip(nor, tum):
        ax.plot([0, 1], [a, b], color=GREY, lw=0.8, alpha=0.55, zorder=1)
    rng = np.random.default_rng(0)
    ax.plot(rng.uniform(-0.045, 0.045, nor.size), nor, "o", color=BLUE, ms=6, zorder=3)
    ax.plot(1 + rng.uniform(-0.045, 0.045, tum.size), tum, "o", color=ORANGE, ms=6, zorder=3)
    ax.plot([0, 1], [nor.mean(), tum.mean()], color="black", lw=3.0, zorder=4)
    ax.set_xlim(-0.28, 1.28)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["adjacent\nnormal", "tumor"])
    ax.set_ylabel("H. pylori log-abundance")
    ax.text(0.42, max(tum.max(), nor.max()) * 0.94,
            f"gastric cancer\n{T.shape[0]} matched pairs\np = {p:.1e}", fontsize=8.5)
    ax.set_title("b  A known carcinogen is DEPLETED in tumor", loc="left", fontsize=10)
    style(ax)
    return diff, p, T.shape[0]


RULES = ["reported", "r>0.99", "r>0.95", "r>0.90", "r>0.80", "ancestor"]
RULE_LABELS = ["reported\n(no collapse)", "r>0.99", "r>0.95", "r>0.90",
               "r>0.80", "ancestor\nrule"]


def figure5(out, dpi):
    """Discovery counts as a function of the taxonomic redundancy rule."""
    import wp0_nesting as nst
    X, meta = nst.load()
    tax = nst.tax_table()

    series = {}
    for label, proj, paired, disp in [("STAD paired", "STAD", True, "STAD (gastric)"),
                                      ("HNSC unpaired", "HNSC", False, "HNSC (head & neck)")]:
        hids, M, allids = nsw.sig_set(X, meta, proj, paired)
        Mh = M[:, [allids.index(h) for h in hids]]
        vals = [len(hids)]
        for t in nsw.R_THRESHOLDS:
            g, _ = nst.redundancy_graph(hids, tax, Mh, rthresh=t)
            vals.append(len(g))
        vals.append(len(nsw.ancestor_clades(hids, tax)))
        series[disp] = vals
        for k, v in zip(RULES, vals):
            exp = nsw.EXPECTED_DISCOVERY[label][k]
            if v != exp:
                DRIFT.append(f"DRIFT F5 {label} {k}: got {v}, expected {exp}")

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    fig.suptitle("Figure 5  Discovery counts depend on the redundancy rule",
                 fontsize=11, x=0.02, ha="left", y=0.985)
    x = range(len(RULES))
    for (name, vals), c, mk in zip(series.items(), [ORANGE, BLUE], ["s", "o"]):
        ax.plot(x, vals, marker=mk, color=c, lw=2.0, ms=7, label=name)
        ax.text(0.06, vals[0], str(vals[0]), fontsize=8.5, color=c, va="center")
        ax.text(len(RULES) - 1.06, vals[-1], str(vals[-1]), fontsize=8.5,
                color=c, va="center", ha="right")

    hn = series["HNSC (head & neck)"]
    fold = hn[0] / max(hn[-1], 1)
    ax.annotate(f"{fold:.1f}-fold\nsame data", xy=(len(RULES) - 1, hn[-1]),
                xytext=(3.05, hn[0] * 0.62), fontsize=9, color=ORANGE,
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))

    ax.set_xticks(list(x))
    ax.set_xticklabels(RULE_LABELS, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel("number of significant findings")
    ax.set_ylim(0, max(hn) * 1.12)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    style(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(out, "Figure5_nesting.png"), dpi=dpi)
    plt.close(fig)
    return series


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()
    if not os.environ.get("TCMA_DIR"):
        print("TCMA_DIR not set. See docs/RUNBOOK.md for the download.", file=sys.stderr)
        return 2
    os.makedirs(a.outdir, exist_ok=True)

    X, meta = tvn.load()
    fig, ax = plt.subplots(1, 2, figsize=(9.55, 5.0))
    fig.suptitle("Figure 4  In real decontaminated data, signal is tissue-of-origin "
                 "and the one causal organism runs backwards",
                 fontsize=11, x=0.01, ha="left", y=0.985)
    aucs, ns = panel_a(ax[0], X, meta)
    diff, p, npairs = panel_b(ax[1], X, meta)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(os.path.join(a.outdir, "Figure4_real_data.png"), dpi=a.dpi)
    plt.close(fig)

    print("AUC: " + "  ".join(f"{p_}={v:.3f} (n={t}v{n})"
                              for p_, v, (t, n) in zip(PROJECTS, aucs, ns)))
    print(f"H. pylori paired diff {diff:+.4f}, p = {p:.2e}, {npairs} pairs")

    series = figure5(a.outdir, a.dpi)
    for name, vals in series.items():
        print(f"Figure 5 {name}: " + " -> ".join(str(v) for v in vals))
    if DRIFT:
        print("\n".join(DRIFT), file=sys.stderr)
        return 1
    print("All asserted values re-derived within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
