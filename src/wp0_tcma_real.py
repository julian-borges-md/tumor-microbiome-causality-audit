"""
RO-2026-008 | WP0 | SIA_001
Analysis 4: FIRST APPLICATION TO REAL DATA.

Dataset: The Cancer Microbiome Atlas (TCMA), Dohlman et al. Cell Host & Microbe
2021, DOI 10.7924/r4bk1j35s. WGS, tissue, sample-level, decontaminated.
This is the field's reference DECONTAMINATED resource, i.e. the best case.

Test: does cancer-type-discriminative signal in TCMA survive the audit that
zero-signal synthetic data failed?
  T1  beats no-information rate
  T3  beats batch-only (sequencing center) baseline
  T5a survives within-batch cross-validation
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import OneHotEncoder

BASE = "/home/claude/tcma"
N_TREES = 150


def load():
    X = pd.read_csv(f"{BASE}/wgs/sample/bacteria.unambiguous.decontam.tissue.sample.rpm.txt",
                    sep="\t", index_col=0, low_memory=False)
    X = X.T                                     # samples x taxa
    sm = pd.read_csv(f"{BASE}/meta/metadata.TCMA.sample.txt", sep="\t", low_memory=False)
    fm = pd.read_csv(f"{BASE}/meta/metadata.TCMA.file.txt", sep="\t", low_memory=False)

    sm = sm[["bcr_sample_barcode", "project", "sample_type"]].dropna(subset=["bcr_sample_barcode"])
    sm = sm.drop_duplicates("bcr_sample_barcode").set_index("bcr_sample_barcode")

    ctr_col = "sequencing_center"
    bcol = "aliquot.bcr_sample_barcode" if "aliquot.bcr_sample_barcode" in fm.columns else None
    fb = fm[[bcol, ctr_col]].dropna().drop_duplicates(bcol).set_index(bcol)

    df = pd.DataFrame(index=X.index)
    df["project"] = sm["project"].reindex(df.index)
    df["center"] = fb[ctr_col].reindex(df.index)
    df["sample_type"] = sm["sample_type"].reindex(df.index)
    keep = df["project"].notna() & df["center"].notna()
    return X.loc[keep], df.loc[keep]


def audit(X, y, batch, label):
    y = pd.factorize(y)[0]
    b = pd.factorize(batch)[0]
    Xv = np.log1p(X.values.astype(float))
    Xv = Xv[:, Xv.sum(0) > 0]
    nir = np.bincount(y).max() / y.size
    cv = StratifiedKFold(5, shuffle=True, random_state=0)

    def clf():
        return RandomForestClassifier(n_estimators=N_TREES, n_jobs=-1,
                                      random_state=0, min_samples_leaf=2)

    acc = cross_val_score(clf(), Xv, y, cv=cv, n_jobs=1).mean()
    p = stats.binomtest(int(round(acc * y.size)), y.size, nir, alternative="greater").pvalue
    t1 = (acc - nir > 0.05) and (p < 0.05)

    enc = OneHotEncoder(sparse_output=False).fit_transform(b.reshape(-1, 1))
    bacc = cross_val_score(clf(), enc, y, cv=cv, n_jobs=1).mean()
    t3 = (acc - bacc) > 0.10

    wb, wn, det = [], [], []
    for bb in np.unique(b):
        m = b == bb
        if m.sum() < 40:
            continue
        yb = y[m]
        keep = np.isin(yb, [c for c in np.unique(yb) if (yb == c).sum() >= 10])
        if len(np.unique(yb[keep])) < 2 or keep.sum() < 40:
            continue
        Xb, ybb = Xv[m][keep], yb[keep]
        a = cross_val_score(clf(), Xb, ybb,
                            cv=StratifiedKFold(3, shuffle=True, random_state=0),
                            n_jobs=1).mean()
        n_ = np.bincount(ybb).max() / ybb.size
        wb.append(a); wn.append(n_); det.append((bb, keep.sum(), len(np.unique(ybb)), a, n_))
    if wb:
        wba, wbn = float(np.mean(wb)), float(np.mean(wn))
        t5a = (wba - wbn) > 0.05
    else:
        wba = wbn = float("nan"); t5a = None

    print(f"\n===== {label} =====")
    print(f"  n samples={y.size}  n taxa={Xv.shape[1]}  classes={len(np.unique(y))}  centers={len(np.unique(b))}")
    print(f"  Random-CV accuracy      : {acc:.3f}   NIR={nir:.3f}   delta={acc-nir:+.3f}")
    print(f"  T1 no-information rate  : {'PASS' if t1 else 'FAIL'}  (p={p:.2e})")
    print(f"  T3 batch-only baseline  : {bacc:.3f}  margin={acc-bacc:+.3f}  -> {'PASS' if t3 else 'FAIL'}")
    if wb:
        print(f"  T5a within-batch CV     : {wba:.3f}  wbNIR={wbn:.3f}  delta={wba-wbn:+.3f} -> {'PASS' if t5a else 'FAIL'}")
        for d in det:
            print(f"        center {str(d[0])[:8]}  n={d[1]:4d} classes={d[2]}  acc={d[3]:.3f} nir={d[4]:.3f} delta={d[3]-d[4]:+.3f}")
    else:
        print("  T5a within-batch CV     : NOT EVALUABLE (no center has >=2 adequately sized classes)")
    verdict = t1 and t3 and (t5a is True)
    print(f"  AUDIT VERDICT: {'PASS' if verdict else 'FAIL / INDETERMINATE'}")
    return dict(acc=acc, nir=nir, bacc=bacc, wba=wba, wbn=wbn, t1=t1, t3=t3, t5a=t5a)


if __name__ == "__main__":
    X, meta = load()
    print("loaded:", X.shape, "| projects:", meta['project'].value_counts().to_dict())
    print("| centers:", meta['center'].nunique(), "| sample types:", meta['sample_type'].value_counts().to_dict())
    audit(X, meta["project"], meta["center"], "TCMA WGS tissue, cancer-type discrimination")
