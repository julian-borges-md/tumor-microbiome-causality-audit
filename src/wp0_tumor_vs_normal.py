"""
RO-2026-008 | WP0 | SIA_001
Analysis 5: tumor vs adjacent-normal WITHIN cancer type (TCMA, real data).

Analysis 4 showed cancer-type discrimination passes the audit, but five
cancer types span five different organs, so much of that signal is simply
tissue-of-origin flora (oral vs gastric vs colonic). That is real but not
causally informative.

The causally relevant question: within one tissue, does the TUMOR carry a
different microbial composition than patient-matched adjacent NORMAL tissue?
Tissue of origin is held constant, so anatomy cannot explain a difference.

Reports per cancer type: audit-style classification (tumor vs normal) plus
per-taxon differential abundance with BH-FDR control.
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score

BASE = "/home/claude/tcma"
N_TREES = 300


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
    keep = d["project"].notna() & d["sample_type"].notna()
    return X.loc[keep], d.loc[keep]


def bh(p, q=0.05):
    p = np.asarray(p); o = np.argsort(p); m = p.size
    thr = q * np.arange(1, m + 1) / m
    ok = p[o] <= thr
    out = np.zeros(m, bool)
    if ok.any():
        out[o[:np.max(np.where(ok)[0]) + 1]] = True
    return out


def taxnames():
    t = pd.read_csv(f"{BASE}/meta/taxonomy.txt", sep="\t", low_memory=False)
    t["tax_id"] = t["tax_id"].astype(float)
    d = dict(zip(t["tax_id"], t["name"]))
    return {**d, **{str(int(k)): v for k, v in d.items() if pd.notna(k)}}


if __name__ == "__main__":
    X, meta = load()
    tn = taxnames()
    print("cohort:", X.shape)
    ct = pd.crosstab(meta["project"], meta["sample_type"])
    print(ct, "\n")

    for proj in ["COAD", "STAD", "HNSC", "ESCA", "READ"]:
        m = meta["project"] == proj
        sub, subm = X.loc[m], meta.loc[m]
        st = subm["sample_type"]
        tum = st.str.contains("Tumor", na=False)
        nor = st.str.contains("Normal", na=False)
        if tum.sum() < 15 or nor.sum() < 10:
            print(f"--- {proj}: skipped (tumor={tum.sum()}, normal={nor.sum()}) ---")
            continue
        sel = tum | nor
        Xv = np.log1p(sub.loc[sel].values.astype(float))
        keepc = Xv.sum(0) > 0
        Xv = Xv[:, keepc]
        cols = sub.columns[keepc]
        y = tum[sel].astype(int).values
        nir = max(y.mean(), 1 - y.mean())

        acc = cross_val_score(
            RandomForestClassifier(n_estimators=N_TREES, n_jobs=-1,
                                   random_state=0, min_samples_leaf=2),
            Xv, y, cv=StratifiedKFold(5, shuffle=True, random_state=0), n_jobs=1).mean()
        p_acc = stats.binomtest(int(round(acc * y.size)), y.size, nir,
                                alternative="greater").pvalue

        pv = np.array([stats.mannwhitneyu(Xv[y == 1, j], Xv[y == 0, j],
                                          alternative="two-sided").pvalue
                       for j in range(Xv.shape[1])])
        pv = np.nan_to_num(pv, nan=1.0)
        hit = bh(pv, 0.05)

        print(f"--- {proj}: tumor={int(y.sum())} normal={int((1-y).sum())} taxa={Xv.shape[1]} ---")
        print(f"    classification acc={acc:.3f}  NIR={nir:.3f}  delta={acc-nir:+.3f}  p={p_acc:.2e}"
              f"   -> {'SIGNAL' if (acc-nir>0.05 and p_acc<0.05) else 'no signal above threshold'}")
        print(f"    differentially abundant taxa at BH-FDR<0.05: {hit.sum()} / {Xv.shape[1]}")
        if hit.sum():
            ordr = np.argsort(pv)
            shown = 0
            for j in ordr:
                if not hit[j]:
                    continue
                lfc = Xv[y == 1, j].mean() - Xv[y == 0, j].mean()
                raw = cols[j]
                try:
                    nm = tn.get(float(raw), tn.get(str(raw), str(raw)))
                except Exception:
                    nm = str(raw)
                print(f"        {str(nm)[:46]:46s} log-diff={lfc:+.3f}  p={pv[j]:.2e}")
                shown += 1
                if shown >= 8:
                    break
        print()
