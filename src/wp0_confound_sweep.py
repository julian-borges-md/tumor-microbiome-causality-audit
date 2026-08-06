"""
RO-2026-008 | WP0 | SIA_001
Analysis 2: Confounding sensitivity sweep.

Question: how much apparent predictive signal does batch-outcome confounding
manufacture from data containing ZERO biological signal, and at what
confounding strength does a conventional analysis become uninterpretable?

Sweep CONFOUND_STRENGTH from chance (0.333, no confounding) to near-perfect
(0.95). At each level, run both cohorts and record:
  - random-CV accuracy (what a conventional paper would report)
  - batch-only baseline (T3)
  - within-batch CV (T5a, the corrected decisive test)
"""
import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import OneHotEncoder

N_CANCERS, N_BATCHES = 6, 3
N_TAXA, N_PER_CANCER = 150, 80
N_SIGNAL, EFFECT, N_TREES = 20, 1.20, 50


def simulate(scenario, confound, rng):
    y = np.repeat(np.arange(N_CANCERS), N_PER_CANCER)
    n = y.size
    home = y // 2
    batch = np.empty(n, int)
    for i in range(n):
        batch[i] = home[i] if rng.random() < confound else \
            rng.choice([b for b in range(N_BATCHES) if b != home[i]])
    base = rng.normal(3.0, 1.0, (n, N_TAXA))
    sig = np.arange(N_SIGNAL)
    bat = np.arange(N_SIGNAL, 2 * N_SIGNAL)
    if scenario == "A":
        for k, t in enumerate(sig):
            base[y == (k % N_CANCERS), t] += EFFECT
    for k, t in enumerate(bat):
        base[batch == (k % N_BATCHES), t] += EFFECT
    return np.log1p(rng.poisson(np.exp(base))), y, batch


def clf():
    return RandomForestClassifier(n_estimators=N_TREES, n_jobs=-1,
                                  random_state=0, min_samples_leaf=2)


def metrics(scenario, confound):
    rng = np.random.default_rng(20260722)
    X, y, batch = simulate(scenario, confound, rng)
    cv = StratifiedKFold(3, shuffle=True, random_state=0)
    nir = np.bincount(y).max() / y.size
    acc = cross_val_score(clf(), X, y, cv=cv, n_jobs=1).mean()
    enc = OneHotEncoder(sparse_output=False).fit_transform(batch.reshape(-1, 1))
    b_acc = cross_val_score(clf(), enc, y, cv=cv, n_jobs=1).mean()
    wb, wn = [], []
    for b in range(N_BATCHES):
        m = batch == b
        yb = y[m]
        keep = np.isin(yb, [c for c in np.unique(yb) if (yb == c).sum() >= 8])
        if len(np.unique(yb[keep])) < 2:
            continue
        Xb, ybb = X[m][keep], yb[keep]
        wb.append(cross_val_score(clf(), Xb, ybb,
                  cv=StratifiedKFold(3, shuffle=True, random_state=0),
                  n_jobs=1).mean())
        wn.append(np.bincount(ybb).max() / ybb.size)
    return dict(conf=confound, scen=scenario, nir=nir, acc=acc,
                batch_only=b_acc, wb=float(np.mean(wb)), wbnir=float(np.mean(wn)))


if __name__ == "__main__":
    scen = sys.argv[1]
    print(f"scenario={scen}  ('A'=real signal, 'B'=ZERO biological signal)")
    print(f"{'confound':>9} {'randomCV':>9} {'NIR':>7} {'inflation':>10} "
          f"{'batchonly':>10} {'withinCV':>9} {'wbNIR':>7} {'wbDelta':>8}")
    for c in [0.333, 0.50, 0.65, 0.80, 0.95]:
        r = metrics(scen, c)
        print(f"{r['conf']:>9.3f} {r['acc']:>9.3f} {r['nir']:>7.3f} "
              f"{r['acc']-r['nir']:>+10.3f} {r['batch_only']:>10.3f} "
              f"{r['wb']:>9.3f} {r['wbnir']:>7.3f} {r['wb']-r['wbnir']:>+8.3f}")
