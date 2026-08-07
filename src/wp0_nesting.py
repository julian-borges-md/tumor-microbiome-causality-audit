"""
RO-2026-008 | WP0 | SIA_001
Analysis 7: taxonomic nesting redundancy.

Motivation. In Analysis 6, 8 "significant" STAD taxa collapsed to 2 organisms
once nested ranks (species / genus / family / order / class) were merged.
Taxonomic hierarchies produce STRUCTURALLY collinear features: a genus column
is, by construction, the sum of its species columns. This has two consequences
nobody controls for consistently:

  (1) DISCOVERY INFLATION: one biological finding is reported as k findings,
      where k = number of ranks at which it reaches significance.
  (2) FDR DISTORTION: BH assumes independence or positive regression
      dependence across tests. Perfectly nested features violate the spirit
      of the correction: the nominal test count overstates the number of
      independent hypotheses actually being asked.

This script measures both, on real TCMA data.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

# TCMA download location. Override with the TCMA_DIR environment
# variable; defaults to ./tcma relative to the working directory.
BASE = os.environ.get("TCMA_DIR", os.path.abspath("tcma"))


def load():
    X = pd.read_csv(f"{BASE}/wgs/sample/bacteria.unambiguous.decontam.tissue.sample.rpm.txt",
                    sep="\t", index_col=0, low_memory=False).T
    sm = pd.read_csv(f"{BASE}/meta/metadata.TCMA.sample.txt", sep="\t", low_memory=False)
    sm = sm[["bcr_sample_barcode", "project", "sample_type", "bcr_patient_barcode"]]
    sm = sm.dropna(subset=["bcr_sample_barcode"]).drop_duplicates("bcr_sample_barcode")
    sm = sm.set_index("bcr_sample_barcode")
    d = pd.DataFrame(index=X.index)
    for c in ["project", "sample_type", "bcr_patient_barcode"]:
        d[c] = sm[c].reindex(d.index)
    keep = d["project"].notna() & d["sample_type"].notna() & d["bcr_patient_barcode"].notna()
    return X.loc[keep], d.loc[keep]


def tax_table():
    t = pd.read_csv(f"{BASE}/meta/taxonomy.txt", sep="\t", low_memory=False)
    t["tax_id"] = t["tax_id"].astype(float)
    t["lineage"] = t["taxonomy"].astype(str)
    return t.set_index("tax_id")


def bh(p, q=0.05):
    p = np.asarray(p); o = np.argsort(p); m = p.size
    thr = q * np.arange(1, m + 1) / m
    ok = p[o] <= thr
    out = np.zeros(m, bool)
    if ok.any():
        out[o[:np.max(np.where(ok)[0]) + 1]] = True
    return out


def redundancy_graph(ids, tax, M, rthresh=0.99):
    """Ancestor-descendant pairs that are near-perfectly collinear.
    A higher rank is redundant with a lower one only if it adds no
    information (|r| > rthresh), not merely because it is an ancestor."""
    lin = {}
    for i in ids:
        try:
            lin[i] = tax.loc[float(i), "lineage"]
        except Exception:
            lin[i] = f"UNKNOWN_{i}"
    idx = {v: k for k, v in enumerate(ids)}
    parent = {i: i for i in ids}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    order = sorted(ids, key=lambda z: len(lin[z]))
    npairs = 0
    for a in range(len(order)):
        la = lin[order[a]]
        for b in range(a + 1, len(order)):
            lb = lin[order[b]]
            if not lb.startswith(la):
                continue
            va, vb = M[:, idx[order[a]]], M[:, idx[order[b]]]
            if va.std() > 0 and vb.std() > 0:
                if abs(np.corrcoef(va, vb)[0, 1]) > rthresh:
                    union(order[a], order[b]); npairs += 1
    groups = {}
    for i in ids:
        groups.setdefault(find(i), []).append(i)
    return list(groups.values()), npairs


def collapse_nested(ids, tax, M, rthresh=0.99):
    g, _ = redundancy_graph(ids, tax, M, rthresh)
    return g


def paired_matrix(X, meta, proj):
    m = meta["project"] == proj
    sub, subm = X.loc[m], meta.loc[m]
    tum = subm["sample_type"].str.contains("Tumor", na=False)
    nor = subm["sample_type"].str.contains("Normal", na=False)
    pt, pn = {}, {}
    for pid, g in subm.groupby("bcr_patient_barcode"):
        ti = g.index[tum.reindex(g.index, fill_value=False)]
        ni = g.index[nor.reindex(g.index, fill_value=False)]
        if len(ti) and len(ni):
            pt[pid] = np.log1p(sub.loc[ti].values.astype(float)).mean(0)
            pn[pid] = np.log1p(sub.loc[ni].values.astype(float)).mean(0)
    if not pt:
        return None, None, None
    pids = sorted(pt)
    return (np.vstack([pt[p] for p in pids]), np.vstack([pn[p] for p in pids]),
            np.array(sub.columns))


if __name__ == "__main__":
    X, meta = load()
    tax = tax_table()

    print("=== PART A: structural redundancy of the feature space ===\n")
    print(f"{'cancer':>7} {'nominal':>8} {'clades':>7} {'redundancy':>11} {'r>0.99 pairs':>13}")
    for proj in ["COAD", "STAD", "HNSC", "ESCA"]:
        T, N, cols = paired_matrix(X, meta, proj)
        if T is None:
            continue
        nz = (T.sum(0) + N.sum(0)) > 0
        ids = cols[nz]
        M = np.vstack([T[:, nz], N[:, nz]])
        clades, npc = redundancy_graph(list(ids), tax, M)
        print(f"{proj:>7} {len(ids):>8} {len(clades):>7} {len(ids)/max(len(clades),1):>10.2f}x {npc:>13}")

    print("\n=== PART B: discovery inflation in the actual results ===\n")
    for proj in ["COAD", "STAD", "HNSC", "ESCA"]:
        T, N, cols = paired_matrix(X, meta, proj)
        if T is None or T.shape[0] < 5:
            continue
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
        pv = np.nan_to_num(pv, nan=1.0)
        hit = bh(pv, 0.05)
        hids = list(ids[hit])
        if not hids:
            print(f"{proj:>7}: 0 reported -> 0 independent  (inflation n/a)")
            continue
        Mh = np.vstack([Ts, Ns])
        cl = collapse_nested(hids, tax, Mh)
        names = []
        for c in cl:
            try:
                names.append(str(tax.loc[float(c[0]), "name"]))
            except Exception:
                names.append(str(c[0]))
        print(f"{proj:>7}: {len(hids)} reported -> {len(cl)} independent clades "
              f"(inflation {len(hids)/len(cl):.1f}x)")
        for c, nm in zip(cl, names):
            print(f"           clade of {len(c):2d} rank(s), representative: {nm}")

    print("\n=== PART C: unpaired HNSC result from Analysis 5 (51 reported hits) ===")
    proj = "HNSC"
    m = meta["project"] == proj
    sub, subm = X.loc[m], meta.loc[m]
    tum = subm["sample_type"].str.contains("Tumor", na=False)
    nor = subm["sample_type"].str.contains("Normal", na=False)
    sel = tum | nor
    Xv = np.log1p(sub.loc[sel].values.astype(float))
    keepc = Xv.sum(0) > 0
    Xv, ids = Xv[:, keepc], np.array(sub.columns)[keepc]
    y = tum[sel].astype(int).values
    pv = np.array([stats.mannwhitneyu(Xv[y == 1, j], Xv[y == 0, j],
                                      alternative="two-sided").pvalue
                   for j in range(Xv.shape[1])])
    pv = np.nan_to_num(pv, nan=1.0)
    hit = bh(pv, 0.05)
    hids = list(ids[hit])
    cl = collapse_nested(hids, tax, Xv)
    print(f"   {len(hids)} reported -> {len(cl)} independent clades "
          f"(inflation {len(hids)/max(len(cl),1):.1f}x)")
    for c in cl:
        try:
            nm = str(tax.loc[float(c[0]), "name"])
        except Exception:
            nm = str(c[0])
        print(f"      clade of {len(c):2d}: {nm}")
