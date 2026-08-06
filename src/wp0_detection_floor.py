"""
RO-2026-008 | WP0 | SIA_001
Analysis 3: Detectability floor (statistical limit of detection).

Question: what is the smallest true taxon-cancer effect this pipeline can
recover at a given sample size? Without this, a low signal-survival rate is
uninterpretable, because "nothing survived" and "we could not have seen it"
look identical.

Design: scenario A (real signal) at realistic confounding (0.75), sweeping
true effect size. 20 of 150 taxa carry genuine cancer-linked signal; 130 are
null. Report sensitivity (true signal taxa recovered at BH-FDR<0.05), the
empirical false discovery among null taxa, and within-batch CV delta.
"""
import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score

N_CANCERS, N_BATCHES = 6, 3
N_TAXA, N_PER_CANCER = 150, 80
N_SIGNAL, CONFOUND, N_TREES = 20, 0.75, 50


def bh(p, q=0.05):
    p = np.asarray(p)
    o = np.argsort(p)
    m = p.size
    thr = q * (np.arange(1, m + 1)) / m
    below = p[o] <= thr
    if not below.any():
        return np.zeros(m, bool)
    kmax = np.max(np.where(below)[0])
    out = np.zeros(m, bool)
    out[o[:kmax + 1]] = True
    return out


def run(effect):
    rng = np.random.default_rng(20260722)
    y = np.repeat(np.arange(N_CANCERS), N_PER_CANCER)
    n = y.size
    home = y // 2
    batch = np.array([home[i] if rng.random() < CONFOUND else
                      rng.choice([b for b in range(N_BATCHES) if b != home[i]])
                      for i in range(n)])
    base = rng.normal(3.0, 1.0, (n, N_TAXA))
    sig = np.arange(N_SIGNAL)
    bat = np.arange(N_SIGNAL, 2 * N_SIGNAL)
    for k, t in enumerate(sig):
        base[y == (k % N_CANCERS), t] += effect
    for k, t in enumerate(bat):
        base[batch == (k % N_BATCHES), t] += 1.20   # batch effect held constant
    X = np.log1p(rng.poisson(np.exp(base)))

    # differential abundance across cancer types, per taxon
    p = np.array([stats.kruskal(*[X[y == c, t] for c in range(N_CANCERS)]).pvalue
                  for t in range(N_TAXA)])
    hit = bh(p, 0.05)
    sens = hit[sig].mean()
    null_idx = np.setdiff1d(np.arange(N_TAXA), np.concatenate([sig, bat]))
    fdr_emp = hit[null_idx].mean()

    # within-batch CV delta
    wb, wn = [], []
    for b in range(N_BATCHES):
        m = batch == b
        yb = y[m]
        keep = np.isin(yb, [c for c in np.unique(yb) if (yb == c).sum() >= 8])
        if len(np.unique(yb[keep])) < 2:
            continue
        Xb, ybb = X[m][keep], yb[keep]
        wb.append(cross_val_score(
            RandomForestClassifier(n_estimators=N_TREES, n_jobs=-1,
                                   random_state=0, min_samples_leaf=2),
            Xb, ybb, cv=StratifiedKFold(3, shuffle=True, random_state=0),
            n_jobs=1).mean())
        wn.append(np.bincount(ybb).max() / ybb.size)
    return sens, fdr_emp, float(np.mean(wb)) - float(np.mean(wn))


if __name__ == "__main__":
    print(f"n={N_CANCERS*N_PER_CANCER} samples, {N_SIGNAL}/{N_TAXA} taxa carry true signal, "
          f"confounding={CONFOUND}")
    print(f"{'effect':>7} {'sensitivity':>12} {'falseDisc':>10} {'wbDelta':>9} {'verdict':>12}")
    for e in [0.10, 0.20, 0.40, 0.60, 0.80, 1.20]:
        s, f, d = run(e)
        v = "detectable" if s >= 0.80 else ("partial" if s >= 0.30 else "BELOW FLOOR")
        print(f"{e:>7.2f} {s:>12.2f} {f:>10.3f} {d:>+9.3f} {v:>12}")
