"""
RO-2026-008 | WP0 | SIA_001
wp0_tcma.py -- canonical real-data analysis with VERIFICATION assertions.

Reproduces every real-data number reported in the session and asserts it
against the recorded value. Any drift fails loudly rather than silently.

Data: The Cancer Microbiome Atlas, DOI 10.7924/r4bk1j35s
      SHA256 of WGS.zip: ec038a07b3b910caa31df07a8f96453f9f6d1a19503504cc43dbc0a76753a333
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import OneHotEncoder

BASE = "/home/claude/tcma"
TREES = 150
TOL = 0.02          # absolute tolerance on accuracies
RPT = {             # values reported during the session
    "n_samples_audit": 611,
    "audit_acc": 0.840,
    "audit_nir": 0.290,
    "audit_batch_only": 0.617,
    "audit_wb_acc": 0.904,
    "stad_pairs": 39,
    "coad_pairs": 21,
    "hnsc_pairs": 22,
    "esca_pairs": 22,
    "hpylori_diff": -0.994,
    "hpylori_p": 4.46e-04,
}
FAILS = []


def check(name, got, exp, tol):
    ok = abs(got - exp) <= tol
    print(f"  {'OK ' if ok else 'DRIFT'}  {name:24s} got={got:<12.5g} reported={exp:<12.5g}")
    if not ok:
        FAILS.append((name, got, exp))
    return ok


def load():
    X = pd.read_csv(f"{BASE}/wgs/sample/bacteria.unambiguous.decontam.tissue.sample.rpm.txt",
                    sep="\t", index_col=0, low_memory=False).T
    sm = pd.read_csv(f"{BASE}/meta/metadata.TCMA.sample.txt", sep="\t", low_memory=False)
    fm = pd.read_csv(f"{BASE}/meta/metadata.TCMA.file.txt", sep="\t", low_memory=False)
    sm = sm[["bcr_sample_barcode", "project", "sample_type", "bcr_patient_barcode"]]
    sm = sm.dropna(subset=["bcr_sample_barcode"]).drop_duplicates("bcr_sample_barcode")
    sm = sm.set_index("bcr_sample_barcode")
    fb = fm[["aliquot.bcr_sample_barcode", "sequencing_center"]].dropna()
    fb = fb.drop_duplicates("aliquot.bcr_sample_barcode").set_index("aliquot.bcr_sample_barcode")
    d = pd.DataFrame(index=X.index)
    for c in ["project", "sample_type", "bcr_patient_barcode"]:
        d[c] = sm[c].reindex(d.index)
    d["center"] = fb["sequencing_center"].reindex(d.index)
    return X, d


def clf():
    return RandomForestClassifier(n_estimators=TREES, n_jobs=-1,
                                  random_state=0, min_samples_leaf=2)


def analysis4(X, meta):
    print("\n[Analysis 4] cancer-type audit")
    keep = meta["project"].notna() & meta["center"].notna()
    Xs, ms = X.loc[keep], meta.loc[keep]
    y = pd.factorize(ms["project"])[0]
    b = pd.factorize(ms["center"])[0]
    Xv = np.log1p(Xs.values.astype(float))
    Xv = Xv[:, Xv.sum(0) > 0]
    check("n_samples", float(y.size), RPT["n_samples_audit"], 0)
    nir = np.bincount(y).max() / y.size
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    acc = cross_val_score(clf(), Xv, y, cv=cv, n_jobs=1).mean()
    enc = OneHotEncoder(sparse_output=False).fit_transform(b.reshape(-1, 1))
    bacc = cross_val_score(clf(), enc, y, cv=cv, n_jobs=1).mean()
    wb, wn = [], []
    for bb in np.unique(b):
        m = b == bb
        if m.sum() < 40:
            continue
        yb = y[m]
        ok = [c for c in np.unique(yb) if (yb == c).sum() >= 10]
        kp = np.isin(yb, ok)
        if len(ok) < 2 or kp.sum() < 40:
            continue
        wb.append(cross_val_score(clf(), Xv[m][kp], yb[kp],
                                  cv=StratifiedKFold(3, shuffle=True, random_state=0),
                                  n_jobs=1).mean())
        wn.append(np.bincount(yb[kp]).max() / yb[kp].size)
    check("audit_acc", acc, RPT["audit_acc"], TOL)
    check("audit_nir", nir, RPT["audit_nir"], TOL)
    check("audit_batch_only", bacc, RPT["audit_batch_only"], TOL)
    check("audit_wb_acc", float(np.mean(wb)), RPT["audit_wb_acc"], TOL)


def pairs(X, meta, proj):
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
    ids = sorted(pt)
    return (np.vstack([pt[p] for p in ids]), np.vstack([pn[p] for p in ids]),
            np.array(sub.columns))


def analysis6(X, meta):
    print("\n[Analysis 6] paired within-patient")
    for proj, key in [("COAD", "coad_pairs"), ("STAD", "stad_pairs"),
                      ("HNSC", "hnsc_pairs"), ("ESCA", "esca_pairs")]:
        T, N, cols = pairs(X, meta, proj)
        n = 0 if T is None else T.shape[0]
        check(f"{proj}_pairs", float(n), RPT[key], 0)

    tax = pd.read_csv(f"{BASE}/meta/taxonomy.txt", sep="\t", low_memory=False)
    tax["tax_id"] = tax["tax_id"].astype(float)
    hp = tax.loc[tax["name"] == "Helicobacter pylori", "tax_id"].iloc[0]
    T, N, cols = pairs(X, meta, "STAD")
    j = int(np.where(cols.astype(str) == str(int(hp)))[0][0])
    d = T[:, j] - N[:, j]
    p = stats.wilcoxon(d, zero_method="wilcox").pvalue
    check("hpylori_diff", float(d.mean()), RPT["hpylori_diff"], 0.01)
    check("hpylori_log10p", float(np.log10(p)), float(np.log10(RPT["hpylori_p"])), 0.15)


if __name__ == "__main__":
    X, meta = load()
    print(f"loaded matrix {X.shape}")
    analysis4(X, meta)
    analysis6(X, meta)
    print("\n" + "=" * 58)
    if FAILS:
        print(f"VERIFICATION FAILED: {len(FAILS)} value(s) drifted")
        for f in FAILS:
            print(f"   {f[0]}: got {f[1]}, reported {f[2]}")
        sys.exit(1)
    print("VERIFICATION PASSED: all reported real-data values reproduce")
