"""
RO-2026-008 | WP0 | SIA_001
Baseline Audit Suite v2.0 -- VALIDATED

VALIDATION RESULT (ground-truth harness, seed 20260722)
--------------------------------------------------------------------------
                          COHORT A            COHORT B
                       (real signal)     (pure batch artifact,
                                          ZERO real signal)
  Random-CV accuracy      0.677               0.332      [NIR = 0.167]
  T1 no-information rate  PASS (+0.511)       PASS (+0.165)   <-- FAILS TO CATCH
  T2 label permutation    PASS                PASS            <-- FAILS TO CATCH
  T3 batch-only baseline  PASS (+0.288)       FAIL (-0.058)   <-- CATCHES
  T5a within-batch CV     PASS (+0.265)       FAIL (-0.022)   <-- CATCHES
  VERDICT                 PASS                FAIL
--------------------------------------------------------------------------
Audit is fit for purpose: passes genuine signal, rejects pure artifact.

CRITICAL FINDING: T1 and T2, the two significance tests most commonly
reported in the microbiome literature, BOTH PASSED on data containing zero
biological signal. They are necessary but insufficient. Only the
confounder-baseline comparison (T3) and within-batch validation (T5a)
discriminate.

v1 -> v2 CORRECTION: leave-one-batch-out (LOBO) is INVALID as a gate under
strong batch-outcome confounding, because holding out a batch also removes
most of an outcome class from training. v1 rejected genuine signal for this
reason (false negative). LOBO retained as diagnostic only.
"""

import sys
import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import OneHotEncoder

N_CANCERS = 6
N_BATCHES = 3
N_TAXA = 250
N_PER_CANCER = 110
N_SIGNAL_TAXA = 20
EFFECT = 1.20
CONFOUND_STRENGTH = 0.75
N_PERM = 25
N_TREES = 60


def simulate(scenario, rng):
    y = np.repeat(np.arange(N_CANCERS), N_PER_CANCER)
    n = y.size
    home = y // 2
    batch = np.empty(n, dtype=int)
    for i in range(n):
        if rng.random() < CONFOUND_STRENGTH:
            batch[i] = home[i]
        else:
            batch[i] = rng.choice([b for b in range(N_BATCHES) if b != home[i]])

    base = rng.normal(3.0, 1.0, size=(n, N_TAXA))
    sig = np.arange(N_SIGNAL_TAXA)
    bat = np.arange(N_SIGNAL_TAXA, 2 * N_SIGNAL_TAXA)

    if scenario == "A":
        for k, t in enumerate(sig):
            base[y == (k % N_CANCERS), t] += EFFECT
        for k, t in enumerate(bat):
            base[batch == (k % N_BATCHES), t] += EFFECT
    else:  # B: pure artifact, zero biological signal
        for k, t in enumerate(bat):
            base[batch == (k % N_BATCHES), t] += EFFECT

    X = np.log1p(rng.poisson(np.exp(base)))
    return X, y, batch


def clf():
    return RandomForestClassifier(n_estimators=N_TREES, n_jobs=-1,
                                  random_state=0, min_samples_leaf=2)


def run(scenario):
    rng = np.random.default_rng(20260722)
    X, y, batch = simulate(scenario, rng)
    nir = np.bincount(y).max() / y.size

    cv5 = StratifiedKFold(5, shuffle=True, random_state=0)
    acc = cross_val_score(clf(), X, y, cv=cv5, n_jobs=1).mean()

    # T1 no-information rate
    p_nir = stats.binomtest(int(round(acc * y.size)), y.size, nir,
                            alternative="greater").pvalue
    t1 = acc - nir > 0.05 and p_nir < 0.05

    # T2 label permutation
    cv3 = StratifiedKFold(3, shuffle=True, random_state=0)
    null = np.array([cross_val_score(clf(), X, rng.permutation(y),
                                     cv=cv3, n_jobs=1).mean()
                     for _ in range(N_PERM)])
    p95 = np.percentile(null, 95)
    t2 = acc > p95

    # T3 batch-only predictor
    enc = OneHotEncoder(sparse_output=False).fit_transform(batch.reshape(-1, 1))
    b_acc = cross_val_score(clf(), enc, y, cv=cv5, n_jobs=1).mean()
    t3 = (acc - b_acc) > 0.10

    # T5a WITHIN-BATCH cross-validation  <-- decisive corrected test
    wb = []
    for b in range(N_BATCHES):
        m = batch == b
        if m.sum() < 60:
            continue
        yb = y[m]
        keep = np.isin(yb, [c for c in np.unique(yb) if (yb == c).sum() >= 10])
        if len(np.unique(yb[keep])) < 2:
            continue
        Xb, ybb = X[m][keep], yb[keep]
        nb = np.bincount(ybb).max() / ybb.size
        a = cross_val_score(clf(), Xb, ybb,
                            cv=StratifiedKFold(3, shuffle=True, random_state=0),
                            n_jobs=1).mean()
        wb.append((a, nb))
    wb_acc = float(np.mean([a for a, _ in wb]))
    wb_nir = float(np.mean([n for _, n in wb]))
    t5a = (wb_acc - wb_nir) > 0.05

    # T5b LOBO, diagnostic only
    lo = []
    for b in range(N_BATCHES):
        tr, te = batch != b, batch == b
        if te.sum() == 0:
            continue
        m = clf().fit(X[tr], y[tr])
        lo.append((m.predict(X[te]) == y[te]).mean())
    lobo = float(np.mean(lo))

    verdict = t1 and t2 and t3 and t5a
    return dict(scenario=scenario, acc=acc, nir=nir, t1=t1, null_p95=p95, t2=t2,
                batch_only=b_acc, t3=t3, wb_acc=wb_acc, wb_nir=wb_nir, t5a=t5a,
                lobo=lobo, verdict=verdict)


if __name__ == "__main__":
    sc = sys.argv[1]
    r = run(sc)
    lbl = "real signal present" if sc == "A" else "PURE BATCH ARTIFACT, zero real signal"
    print(f"=== COHORT {sc} ({lbl}) ===")
    print(f"  Random-CV accuracy (conventional)      : {r['acc']:.3f}   [NIR={r['nir']:.3f}]")
    print(f"  T1 no-information rate                 : {'PASS' if r['t1'] else 'FAIL'}  delta={r['acc']-r['nir']:+.3f}")
    print(f"  T2 label permutation (null p95={r['null_p95']:.3f}) : {'PASS' if r['t2'] else 'FAIL'}")
    print(f"  T3 batch-only predictor={r['batch_only']:.3f}       : {'PASS' if r['t3'] else 'FAIL'}  margin={r['acc']-r['batch_only']:+.3f}")
    print(f"  T5a WITHIN-BATCH CV={r['wb_acc']:.3f} [nir={r['wb_nir']:.3f}] : {'PASS' if r['t5a'] else 'FAIL'}  delta={r['wb_acc']-r['wb_nir']:+.3f}")
    print(f"  T5b LOBO (diagnostic only)             : {r['lobo']:.3f}")
    print(f"  AUDIT VERDICT: {'PASS - signal trusted' if r['verdict'] else 'FAIL - signal rejected'}")
