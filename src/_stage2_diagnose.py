#!/usr/bin/env python3
"""Diagnostic for RLO Stage-2 validation. NOT the pipeline. Determines which
selection/estimator variant reproduces the committed tables, to diagnose the
committed pipeline's discrepancies. Does not modify wp1_mr_pipeline.py."""
import csv, sys, os, math
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from wp1_mr_pipeline import load_exposure, ivw, cochran_q
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PAL = {frozenset("AT"), frozenset("CG")}

def select(hits, window, pal_before):
    keep = [h for h in hits if h["p"] < 1e-5]
    if pal_before:
        keep = [h for h in keep if frozenset({h["ref"], h["eff"]}) not in PAL]
    keep.sort(key=lambda h: h["p"])
    chosen = []
    for h in keep:
        if window and any(h["chrom"] == c["chrom"] and abs(h["bp"] - c["bp"]) < window
                          for c in chosen):
            continue
        chosen.append(h)
    return chosen

def load_out(path):
    out = {}
    with open(path) as fh:
        hd = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hd)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try: b = float(f[ix["beta"]]); se = float(f[ix["sebeta"]])
            except ValueError: continue
            rec = dict(ea=f[ix["alt"]].upper(), oa=f[ix["ref"]].upper(), beta=b, se=se)
            for rs in f[ix["rsids"]].split(","):
                if rs.startswith("rs"): out[rs] = rec
    return out

def harm(h, o):
    ea, oa, eo, oo = h["eff"], h["ref"], o["ea"], o["oa"]
    if frozenset({ea, oa}) in PAL: return None
    if ea == eo and oa == oo: bo = o["beta"]
    elif ea == oo and oa == eo: bo = -o["beta"]
    elif ea == COMP.get(eo) and oa == COMP.get(oo): bo = o["beta"]
    elif ea == COMP.get(oo) and oa == COMP.get(eo): bo = -o["beta"]
    else: return None
    return h["beta"], bo, o["se"]

def egger_oriented(bx, by, sy):
    s = np.sign(bx); s[s == 0] = 1
    bxo, byo = bx * s, by * s
    w = 1.0 / sy ** 2
    X = np.column_stack([np.ones_like(bxo), bxo]); cov = np.linalg.inv(X.T @ (w[:, None] * X))
    coef = cov @ (X.T @ (w * byo)); resid = byo - X @ coef; dof = max(len(bxo) - 2, 1)
    seb = np.sqrt(np.diag(cov) * (resid ** 2 * w).sum() / dof)
    return 2 * stats.t.sf(abs(coef[1] / seb[1]), dof), 2 * stats.t.sf(abs(coef[0] / seb[0]), dof)

def build(expo, outc, window, pal_before):
    res = {}
    for tax, hits in expo.items():
        bx, by, sy = [], [], []
        for h in select(hits, window, pal_before):
            o = outc.get(h["rsid"])
            if o is None: continue
            z = harm(h, o)
            if z is None: continue
            bx.append(z[0]); by.append(z[1]); sy.append(z[2])
        if len(bx) < 3: continue
        bx, by, sy = np.array(bx), np.array(by), np.array(sy)
        b, se, p = ivw(bx, by, sy)
        ep, eip = egger_oriented(bx, by, sy)
        q, qp = cochran_q(bx, by, sy, b)
        res[tax] = dict(n=len(bx), p=p, ep=ep, eip=eip, qp=qp)
    return res

def approx(a, b, rel=1e-4):
    try: a, b = float(a), float(b)
    except: return False
    if a == 0 and b == 0: return True
    return abs(a - b) <= 1e-9 + rel * abs(a)

if __name__ == "__main__":
    expo = load_exposure("data/raw/MBG.allHits.p1e4.txt")
    outc = load_out("outcome/colorectal.hits.tsv")
    comm = {r["taxon"]: r for r in csv.DictReader(open("results/MR_results_colorectal.tsv"), delimiter="\t")}
    print(f"{'variant':<34}{'n match':>9}{'p_ivw':>9}{'p_egger':>9}{'q_p':>7}")
    for label, window, pb in [("500kb, pal-before (=pipeline)", 500_000, True),
                              ("1Mb, pal-before", 1_000_000, True),
                              ("500kb, pal-after", 500_000, False),
                              ("1Mb, pal-after (=cloud original?)", 1_000_000, False)]:
        r = build(expo, outc, window, pb)
        nm = pm = em = qm = tot = 0
        for t, c in comm.items():
            if t not in r: continue
            tot += 1
            if r[t]["n"] == int(c["n_snp"]): nm += 1
            if approx(c["p_ivw"], r[t]["p"]): pm += 1
            if approx(c["p_egger"], r[t]["ep"]): em += 1
            if approx(c["q_p"], r[t]["qp"]): qm += 1
        N = len(comm)
        print(f"{label:<34}{nm:>4}/{N:<4}{pm:>5}/{N:<3}{em:>5}/{N:<3}{qm:>3}/{N}")
