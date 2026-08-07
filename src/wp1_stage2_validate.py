#!/usr/bin/env python3
"""
RO-2026-008 | WP1 | Stage 2 reconstruction and validation harness
Serves RLO_MR_OUTCOME-STAGE2-VALIDATION_v1.0.

The committed wp1_mr_pipeline.py `run` subcommand is a stub: it prints the
"not validated" notice and exits 2. Only the estimator functions and the
Stage-1 selection are committed. This harness imports those committed
functions unchanged and supplies the outcome harmonisation + per-taxon loop
that `run` never wired up, so Stage 2 can be executed and compared to the
committed tables WITHOUT modifying the pipeline (RLO Section 12).

Harmonisation follows the documented method: match on rsID, drop palindromic
(already dropped in selection), align effect alleles across strand and coding,
drop allele mismatches and instruments absent from the outcome.

Usage:
    python3 src/wp1_stage2_validate.py \
        --exposure data/raw/MBG.allHits.p1e4.txt \
        --outcome-hits outcome/colorectal.hits.tsv \
        --fmt finngen \
        --committed results/MR_results_colorectal.tsv \
        --out results/MR_rebuilt_colorectal.tsv
"""
import argparse
import csv
import sys
import numpy as np
from scipy import stats

sys.path.insert(0, __import__("os").path.dirname(__file__))
from wp1_mr_pipeline import (               # committed functions, unchanged
    load_exposure, select_instruments, ivw, mr_egger,
    weighted_median, cochran_q, bh, PALINDROMIC, CLUMP_WINDOW)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_outcome_hits(path, fmt):
    """Return {rsid: dict(ea, oa, beta, se, p)} from a pre-extracted subset."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        if fmt == "finngen":
            c_rs, c_ea, c_oa = idx["rsids"], idx["alt"], idx["ref"]
            c_b, c_se = idx["beta"], idx["sebeta"]
            c_p = idx["pval"]
        elif fmt == "gwascat":
            c_rs = idx.get("rsid", idx.get("rs_id"))
            c_ea, c_oa = idx["effect_allele"], idx["other_allele"]
            c_b, c_se = idx["beta"], idx["standard_error"]
            c_p = idx["p_value"]
        else:
            raise SystemExit(f"unknown fmt {fmt}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= c_p:
                continue
            try:
                beta = float(f[c_b]); se = float(f[c_se])
            except ValueError:
                continue
            rec = dict(ea=f[c_ea].upper(), oa=f[c_oa].upper(),
                       beta=beta, se=se, p=f[c_p])
            for rs in f[c_rs].split(","):        # finngen may comma-join rsids
                rs = rs.strip()
                if rs and rs.startswith("rs"):
                    # keep FIRST row at multi-allelic sites; see CORRECTIONS C18
                    out.setdefault(rs, rec)
    return out


def harmonise(instr, orec):
    """Align one exposure instrument to the outcome record.
    Returns (beta_exp, se_exp, beta_out_aligned, se_out) or None to drop."""
    ea, oa = instr["eff"], instr["ref"]           # exposure effect / other
    eo, oo = orec["ea"], orec["oa"]               # outcome effect / other
    if frozenset({ea, oa}) in PALINDROMIC:
        return None
    if ea == eo and oa == oo:
        bo = orec["beta"]
    elif ea == oo and oa == eo:
        bo = -orec["beta"]
    elif ea == COMP.get(eo) and oa == COMP.get(oo):
        bo = orec["beta"]
    elif ea == COMP.get(oo) and oa == COMP.get(eo):
        bo = -orec["beta"]
    else:
        return None
    return instr["beta"], instr["se"], bo, orec["se"]


def ivw_re(bx, by, sy):
    """Random-effects IVW: fixed-effect SE scaled by sqrt(max(1, Q/df))."""
    w = 1.0 / sy ** 2
    b = float(np.sum(w * bx * by) / np.sum(w * bx ** 2))
    se = float(np.sqrt(1.0 / np.sum(w * bx ** 2)))
    q = float(np.sum(w * (by - b * bx) ** 2))
    df = max(len(bx) - 1, 1)
    se *= max(1.0, q / df) ** 0.5
    return b, se, 2 * stats.norm.sf(abs(b / se))


def egger_oriented(bx, by, sy):
    """MR-Egger with exposure effects oriented positive."""
    s = np.sign(bx); s[s == 0] = 1
    bxo, byo = bx * s, by * s
    w = 1.0 / sy ** 2
    X = np.column_stack([np.ones_like(bxo), bxo])
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    coef = cov @ (X.T @ (w * byo))
    resid = byo - X @ coef
    dof = max(len(bxo) - 2, 1)
    seb = np.sqrt(np.diag(cov) * (resid ** 2 * w).sum() / dof)
    return (float(coef[1]), float(seb[1]),
            2 * stats.t.sf(abs(coef[1] / seb[1]), dof),
            2 * stats.t.sf(abs(coef[0] / seb[0]), dof))


def rebuild(exposure_path, hits_path, fmt, window, faithful=False):
    exposure = load_exposure(exposure_path)
    outcome = load_outcome_hits(hits_path, fmt)
    rows = []
    for tax, hits in exposure.items():
        sel = select_instruments(hits, window=window,
                                 drop_palindromic=not faithful)
        bx, by, sy = [], [], []
        for h in sel:
            orec = outcome.get(h["rsid"])
            if orec is None:
                continue
            hz = harmonise(h, orec)
            if hz is None:
                continue
            bx.append(hz[0]); by.append(hz[2]); sy.append(hz[3])
        n = len(bx)
        if n < 3:
            continue
        bx = np.array(bx); by = np.array(by); sy = np.array(sy)
        if faithful:
            b, se, p = ivw_re(bx, by, sy)
            eb, ese, ep, eip = egger_oriented(bx, by, sy)
        else:
            b, se, p = ivw(bx, by, sy)
            eb, ese, ep, _, eip = mr_egger(bx, by, sy)
        wb, wse, wp = weighted_median(bx, by, sy)
        q, qp = cochran_q(bx, by, sy, b)
        rows.append(dict(taxon=tax, n_snp=n, b_ivw=b, se_ivw=se, p_ivw=p,
                         or_ivw=float(np.exp(b)),
                         lci=float(np.exp(b - 1.96 * se)),
                         uci=float(np.exp(b + 1.96 * se)),
                         b_egger=eb, se_egger=ese, p_egger=ep,
                         egger_intercept_p=eip,
                         b_wmedian=wb, se_wmedian=wse, p_wmedian=wp,
                         q_stat=q, q_p=qp))
    # BH-FDR across taxa, on IVW p
    ps = np.array([r["p_ivw"] for r in rows])
    hit = bh(ps)                                  # committed BH (boolean)
    order = np.argsort(ps)
    m = len(ps)
    q_bh = np.empty(m); prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ps[order[i]] * m / (i + 1))
        q_bh[order[i]] = prev
    for r, fq in zip(rows, q_bh):
        r["fdr_ivw"] = float(fq)
    return rows


FIELDS = ["taxon", "n_snp", "b_ivw", "se_ivw", "p_ivw", "or_ivw", "lci", "uci",
          "fdr_ivw", "b_egger", "se_egger", "p_egger", "egger_intercept_p",
          "b_wmedian", "se_wmedian", "p_wmedian", "q_stat", "q_p"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exposure", required=True)
    ap.add_argument("--outcome-hits", required=True)
    ap.add_argument("--fmt", required=True, choices=["finngen", "gwascat"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=None)
    ap.add_argument("--faithful", action="store_true",
                    help="apply the reconstruction that reproduces the committed "
                         "tables: 1Mb window, palindromic dropped at harmonisation, "
                         "random-effects IVW, oriented MR-Egger")
    a = ap.parse_args()
    window = a.window if a.window is not None else (
        1_000_000 if a.faithful else CLUMP_WINDOW)
    rows = rebuild(a.exposure, a.outcome_hits, a.fmt, window, a.faithful)
    rows.sort(key=lambda r: r["p_ivw"])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in FIELDS})
    mode = "faithful reconstruction" if a.faithful else "committed pipeline params"
    print(f"rebuilt {len(rows)} taxa ({mode}, window={window}) -> {a.out}")


if __name__ == "__main__":
    main()
