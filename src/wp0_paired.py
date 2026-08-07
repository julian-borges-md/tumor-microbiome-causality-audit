"""
RO-2026-008 | WP0 | SIA_001
Analysis 6: PAIRED within-patient tumor vs adjacent normal (TCMA, real data).

Analysis 5 was unpaired and underpowered (n_normal = 21 to 39). Inter-patient
variation in microbiome composition is large and swamps the tumor effect.
Pairing each tumor to the SAME patient's adjacent normal removes that
variance entirely and is the largest available power gain from this dataset.

Test: Wilcoxon signed-rank on within-patient paired differences, BH-FDR.
Positive controls re-tested under the paired design.
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


def taxnames():
    t = pd.read_csv(f"{BASE}/meta/taxonomy.txt", sep="\t", low_memory=False)
    t["tax_id"] = t["tax_id"].astype(float)
    return dict(zip(t["tax_id"], t["name"]))


def bh(p, q=0.05):
    p = np.asarray(p); o = np.argsort(p); m = p.size
    thr = q * np.arange(1, m + 1) / m
    ok = p[o] <= thr
    out = np.zeros(m, bool)
    if ok.any():
        out[o[:np.max(np.where(ok)[0]) + 1]] = True
    return out


def build_pairs(X, meta, proj):
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
    return (np.vstack([pt[p] for p in pids]),
            np.vstack([pn[p] for p in pids]),
            list(sub.columns))


if __name__ == "__main__":
    X, meta = load()
    tn = taxnames()
    name2id = {v: k for k, v in tn.items() if isinstance(v, str)}

    print("=== PAIRED WITHIN-PATIENT TUMOR vs ADJACENT NORMAL (TCMA) ===\n")
    summary = []
    for proj in ["COAD", "STAD", "HNSC", "ESCA", "READ"]:
        T, N, cols = build_pairs(X, meta, proj)
        if T is None or T.shape[0] < 5:
            print(f"--- {proj}: {0 if T is None else T.shape[0]} matched pairs, skipped ---\n")
            continue
        npairs = T.shape[0]
        nz = (T.sum(0) + N.sum(0)) > 0
        Ts, Ns = T[:, nz], N[:, nz]
        cs = np.array(cols)[nz]

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
        up = int(((Ts - Ns).mean(0) > 0)[hit].sum())
        dn = int(hit.sum()) - up

        print(f"--- {proj}: {npairs} matched pairs, {Ts.shape[1]} taxa ---")
        print(f"    significant at BH-FDR<0.05: {hit.sum()}   (enriched in tumor: {up}, depleted: {dn})")
        if hit.sum():
            for j in np.argsort(pv)[:8]:
                if not hit[j]:
                    continue
                diff = (Ts[:, j] - Ns[:, j]).mean()
                try:
                    nm = tn.get(float(cs[j]), str(cs[j]))
                except Exception:
                    nm = str(cs[j])
                arrow = "UP in tumor  " if diff > 0 else "DOWN in tumor"
                print(f"        {str(nm)[:42]:42s} {arrow} diff={diff:+.3f} p={pv[j]:.2e}")
        summary.append((proj, npairs, int(hit.sum()), up, dn))

        # positive controls, paired
        pcs = {"COAD": ["Fusobacterium nucleatum", "Fusobacterium", "Bacteroides fragilis"],
               "STAD": ["Helicobacter pylori"],
               "HNSC": ["Fusobacterium nucleatum"],
               "ESCA": ["Fusobacterium nucleatum"]}.get(proj, [])
        for nm in pcs:
            tid = name2id.get(nm)
            if tid is None:
                continue
            idx = np.where(cs.astype(str) == str(int(tid)))[0]
            if not len(idx):
                print(f"      [PC] {nm:30s} absent after zero-filter")
                continue
            j = idx[0]
            d = Ts[:, j] - Ns[:, j]
            if np.count_nonzero(d) < 3:
                print(f"      [PC] {nm:30s} too few nonzero pairs (n={np.count_nonzero(d)})")
                continue
            p = stats.wilcoxon(d, zero_method="wilcox").pvalue
            print(f"      [PC] {nm:30s} paired diff={d.mean():+.3f}  "
                  f"n_pairs_nonzero={np.count_nonzero(d)}  p={p:.3f}")
        print()

    print("=== SUMMARY ===")
    print(f"{'cancer':>7} {'pairs':>6} {'sig':>5} {'up':>4} {'down':>5}")
    for s in summary:
        print(f"{s[0]:>7} {s[1]:>6} {s[2]:>5} {s[3]:>4} {s[4]:>5}")
