#!/usr/bin/env python3
"""RLO_MR_OUTCOME-STAGE2-VALIDATION_v1.0 | per-taxon comparison, all 211 taxa.

Builds each outcome TWO ways and compares both to the committed table:
  A) COMMITTED PIPELINE parameters: 500kb window, palindromic dropped BEFORE
     clumping, fixed-effect IVW SE, MR-Egger without exposure orientation.
  B) FAITHFUL RECONSTRUCTION: 1Mb window, palindromic dropped AFTER clumping
     (during harmonisation), random-effects IVW SE, oriented MR-Egger.

Emits docs/findings/_pertaxon_{outcome}.tsv with one row per taxon.
Diagnostic only. Does not modify wp1_mr_pipeline.py.
"""
import csv, os, sys
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from wp1_mr_pipeline import load_exposure, ivw as ivw_fixed, mr_egger, weighted_median, cochran_q

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


def load_out(path, fmt, keep_first=True):
    out = {}
    with open(path) as fh:
        hd = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        ix = {c: i for i, c in enumerate(hd)}
        if fmt == "finngen":
            crs, cea, coa, cb, cse = ix["rsids"], ix["alt"], ix["ref"], ix["beta"], ix["sebeta"]
        else:
            crs = ix["rsid"] if "rsid" in ix else ix["rs_id"]
            cea, coa = ix["effect_allele"], ix["other_allele"]
            cb, cse = ix["beta"], ix["standard_error"]
        for line in fh:
            f = line.rstrip("\n").rstrip("\r").split("\t")
            if len(f) <= max(crs, cea, coa, cb, cse):
                continue
            try:
                b = float(f[cb]); se = float(f[cse])
            except ValueError:
                continue
            if not np.isfinite(b) or not np.isfinite(se) or se <= 0:
                continue
            rec = dict(ea=f[cea].upper(), oa=f[coa].upper(), beta=b, se=se)
            for rs in f[crs].split(","):
                rs = rs.strip()
                if rs.startswith("rs"):
                    if keep_first:
                        out.setdefault(rs, rec)
                    else:
                        out[rs] = rec
    return out


def harm(h, o):
    ea, oa, eo, oo = h["eff"], h["ref"], o["ea"], o["oa"]
    if frozenset({ea, oa}) in PAL:
        return None
    if ea == eo and oa == oo:
        bo = o["beta"]
    elif ea == oo and oa == eo:
        bo = -o["beta"]
    elif ea == COMP.get(eo) and oa == COMP.get(oo):
        bo = o["beta"]
    elif ea == COMP.get(oo) and oa == COMP.get(eo):
        bo = -o["beta"]
    else:
        return None
    return h["beta"], bo, o["se"]


def ivw_re(bx, by, sy):
    w = 1.0 / sy ** 2
    b = float(np.sum(w * bx * by) / np.sum(w * bx ** 2))
    se = float(np.sqrt(1.0 / np.sum(w * bx ** 2)))
    q = float(np.sum(w * (by - b * bx) ** 2))
    df = max(len(bx) - 1, 1)
    se *= max(1.0, q / df) ** 0.5
    return b, se, 2 * stats.norm.sf(abs(b / se))


def egger_oriented(bx, by, sy):
    s = np.sign(bx); s[s == 0] = 1
    bxo, byo = bx * s, by * s
    w = 1.0 / sy ** 2
    X = np.column_stack([np.ones_like(bxo), bxo])
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    coef = cov @ (X.T @ (w * byo))
    resid = byo - X @ coef
    dof = max(len(bxo) - 2, 1)
    seb = np.sqrt(np.diag(cov) * (resid ** 2 * w).sum() / dof)
    return (float(coef[1]), 2 * stats.t.sf(abs(coef[1] / seb[1]), dof),
            2 * stats.t.sf(abs(coef[0] / seb[0]), dof))


def build(expo, outc, window, pal_before, faithful):
    res = {}
    for tax, hits in expo.items():
        bx, by, sy = [], [], []
        for h in select(hits, window, pal_before):
            o = outc.get(h["rsid"])
            if o is None:
                continue
            z = harm(h, o)
            if z is None:
                continue
            bx.append(z[0]); by.append(z[1]); sy.append(z[2])
        if len(bx) < 3:
            res[tax] = None
            continue
        bx, by, sy = np.array(bx), np.array(by), np.array(sy)
        if faithful:
            b, se, p = ivw_re(bx, by, sy)
            eb, ep, eip = egger_oriented(bx, by, sy)
        else:
            b, se, p = ivw_fixed(bx, by, sy)
            eb, _, ep, _, eip = mr_egger(bx, by, sy)
        wb, wse, wp = weighted_median(bx, by, sy)
        q, qp = cochran_q(bx, by, sy, b)
        res[tax] = dict(n=len(bx), b=b, se=se, p=p, ep=ep, eip=eip,
                        wp=wp, q=q, qp=qp)
    return res


def approx(a, b, rel=1e-4):
    try:
        a = float(a); b = float(b)
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(a) and np.isfinite(b)):
        return None
    if a == 0 and b == 0:
        return True
    return abs(a - b) <= 1e-12 + rel * abs(a)


OUTCOMES = [
    dict(name="colorectal", hits="outcome/colorectal.hits.tsv", fmt="finngen",
         committed="results/MR_results_colorectal.tsv",
         cm=dict(n="n_snp", b=None, se=None, p="p_ivw", ep="p_egger",
                 eip="egger_intercept_p", wp="p_wmedian", qp="q_p")),
    dict(name="colorectal_100k", hits="outcome/colorectal_100k.hits.tsv", fmt="gwascat",
         committed="results/MR_results_colorectal_100k.tsv",
         cm=dict(n="n", b="b", se="se", p="p", ep="p_eg",
                 eip="pleio_p", wp="p_wm", qp="qp")),
    dict(name="pancreatic", hits="outcome/pancreatic.hits.tsv", fmt="finngen",
         committed="results/MR_results_pancreatic.tsv",
         cm=dict(n="n_snp", b="b_ivw", se="se_ivw", p="p_ivw", ep="p_egger",
                 eip="egger_intercept_p", wp="p_wmedian", qp="q_p")),
]


def main():
    expo = load_exposure("data/raw/MBG.allHits.p1e4.txt")
    os.makedirs("docs/findings", exist_ok=True)
    summary = []
    for oc in OUTCOMES:
        outc = load_out(oc["hits"], oc["fmt"])
        comm = {r["taxon"]: r for r in
                csv.DictReader(open(oc["committed"]), delimiter="\t")}
        A = build(expo, outc, 500_000, True, False)     # committed pipeline
        B = build(expo, outc, 1_000_000, False, True)   # faithful reconstruction
        cm = oc["cm"]
        rows, cnt = [], dict(nA=0, pA=0, eA=0, qA=0, nB=0, pB=0, eB=0, qB=0,
                             bB=0, seB=0, tot=0)
        for tax, c in comm.items():
            a, b = A.get(tax), B.get(tax)
            cnt["tot"] += 1
            cn = int(c[cm["n"]])
            cb = c[cm["b"]] if cm["b"] else ""
            cse = c[cm["se"]] if cm["se"] else ""
            r = dict(taxon=tax, committed_n=cn,
                     A_n=(a["n"] if a else ""), B_n=(b["n"] if b else ""),
                     committed_p=c[cm["p"]],
                     A_p=(a["p"] if a else ""), B_p=(b["p"] if b else ""),
                     committed_b=cb, B_b=(b["b"] if b else ""),
                     committed_se=cse, B_se=(b["se"] if b else ""),
                     committed_p_egger=c[cm["ep"]],
                     A_p_egger=(a["ep"] if a else ""), B_p_egger=(b["ep"] if b else ""),
                     committed_q_p=c[cm["qp"]],
                     A_q_p=(a["qp"] if a else ""), B_q_p=(b["qp"] if b else ""))
            if a:
                if a["n"] == cn: cnt["nA"] += 1
                if approx(c[cm["p"]], a["p"]): cnt["pA"] += 1
                if approx(c[cm["ep"]], a["ep"]): cnt["eA"] += 1
                if approx(c[cm["qp"]], a["qp"]): cnt["qA"] += 1
            if b:
                if b["n"] == cn: cnt["nB"] += 1
                if approx(c[cm["p"]], b["p"]): cnt["pB"] += 1
                if approx(c[cm["ep"]], b["ep"]): cnt["eB"] += 1
                if approx(c[cm["qp"]], b["qp"]): cnt["qB"] += 1
                if cb and approx(cb, b["b"]): cnt["bB"] += 1
                if cse and approx(cse, b["se"]): cnt["seB"] += 1
            r["reproduces_A"] = bool(a and a["n"] == cn and approx(c[cm["p"]], a["p"]))
            r["reproduces_B"] = bool(b and b["n"] == cn and approx(c[cm["p"]], b["p"]))
            rows.append(r)
        out = f"docs/findings/_pertaxon_{oc['name']}.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        summary.append((oc["name"], cnt, out))
        T = cnt["tot"]
        print(f"\n=== {oc['name']} (n taxa committed = {T}) -> {out}")
        print(f"  A committed-pipeline params : n_snp {cnt['nA']}/{T}  p_ivw {cnt['pA']}/{T}"
              f"  p_egger {cnt['eA']}/{T}  q_p {cnt['qA']}/{T}")
        print(f"  B faithful reconstruction   : n_snp {cnt['nB']}/{T}  p_ivw {cnt['pB']}/{T}"
              f"  p_egger {cnt['eB']}/{T}  q_p {cnt['qB']}/{T}")
        if cm["b"]:
            print(f"  B full-precision b {cnt['bB']}/{T}   se {cnt['seB']}/{T}")
        else:
            print("  b/se comparison IMPOSSIBLE: committed table stores a rounded "
                  "display string only (CORRECTIONS C7)")
    return summary


if __name__ == "__main__":
    main()
