"""
RO-2026-008 | WP0 | SIA_001
wp0_core.py  --  canonical, reproducible implementation.

Supersedes the exploratory scripts written during development. Those used
THREE different parameter configurations (taxa 500/300/250/150, permutations
200/40/25, trees 200/80/60/50) because they were tuned to fit execution
windows. Reported values from those runs are therefore not mutually
comparable. This module fixes ONE canonical configuration and re-derives
every synthetic result from it.

Design rules:
  - all parameters in CONFIG, versioned
  - seed is an explicit argument, never global state
  - no exec-based imports
  - deterministic given (CONFIG, seed)

Usage:
  python3 wp0_core.py audit   --seeds 0,1,2,3,4
  python3 wp0_core.py sweep   --seeds 0,1,2
  python3 wp0_core.py floor   --seeds 0,1,2
"""
import argparse
import json
import sys

import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import OneHotEncoder

CONFIG = {
    "version": "wp0_core-1.0.0",
    "n_cancers": 6,
    "n_batches": 3,
    "n_taxa": 150,
    "n_per_cancer": 60,
    "n_signal_taxa": 20,
    "effect": 1.20,
    "batch_effect": 1.20,
    "confound_default": 0.75,
    "n_trees": 40,
    "cv_folds": 3,
    "base_log_mean": 3.0,
    "base_log_sd": 1.0,
    "min_class_in_batch": 8,
    "t1_min_delta": 0.05,
    "t3_min_margin": 0.10,
    "t5a_min_delta": 0.05,
    "fdr_q": 0.05,
}


def _rng(seed):
    return np.random.default_rng(1_000_000 + int(seed))


def simulate(scenario, confound, seed, effect=None, cfg=CONFIG):
    """scenario 'A' = real taxon-cancer signal + batch effect
       scenario 'B' = batch effect ONLY, zero biological signal"""
    rng = _rng(seed)
    eff = cfg["effect"] if effect is None else effect
    y = np.repeat(np.arange(cfg["n_cancers"]), cfg["n_per_cancer"])
    n = y.size
    home = y // 2
    u = rng.random(n)
    alt = rng.integers(0, cfg["n_batches"] - 1, size=n)
    batch = np.where(u < confound, home,
                     (home + 1 + alt) % cfg["n_batches"])
    base = rng.normal(cfg["base_log_mean"], cfg["base_log_sd"], (n, cfg["n_taxa"]))
    sig = np.arange(cfg["n_signal_taxa"])
    bat = np.arange(cfg["n_signal_taxa"], 2 * cfg["n_signal_taxa"])
    if scenario == "A":
        for k, t in enumerate(sig):
            base[y == (k % cfg["n_cancers"]), t] += eff
    for k, t in enumerate(bat):
        base[batch == (k % cfg["n_batches"]), t] += cfg["batch_effect"]
    X = np.log1p(rng.poisson(np.exp(base)))
    return X, y, batch


def _clf(cfg=CONFIG):
    return RandomForestClassifier(n_estimators=cfg["n_trees"], n_jobs=-1,
                                  random_state=0, min_samples_leaf=2)


def _cv(cfg=CONFIG):
    return StratifiedKFold(cfg["cv_folds"], shuffle=True, random_state=0)


def audit(X, y, batch, cfg=CONFIG, do_perm=False, n_perm=25, seed=0):
    nir = np.bincount(y).max() / y.size
    acc = cross_val_score(_clf(cfg), X, y, cv=_cv(cfg), n_jobs=1).mean()
    p_nir = stats.binomtest(int(round(acc * y.size)), y.size, nir,
                            alternative="greater").pvalue
    t1 = bool((acc - nir > cfg["t1_min_delta"]) and (p_nir < 0.05))

    enc = OneHotEncoder(sparse_output=False).fit_transform(batch.reshape(-1, 1))
    bacc = cross_val_score(_clf(cfg), enc, y, cv=_cv(cfg), n_jobs=1).mean()
    t3 = bool((acc - bacc) > cfg["t3_min_margin"])

    wb, wn = [], []
    for b in np.unique(batch):
        m = batch == b
        yb = y[m]
        ok = [c for c in np.unique(yb) if (yb == c).sum() >= cfg["min_class_in_batch"]]
        keep = np.isin(yb, ok)
        if len(ok) < 2 or keep.sum() < 20:
            continue
        wb.append(cross_val_score(_clf(cfg), X[m][keep], yb[keep],
                                  cv=_cv(cfg), n_jobs=1).mean())
        wn.append(np.bincount(yb[keep]).max() / yb[keep].size)
    if wb:
        wba, wbn = float(np.mean(wb)), float(np.mean(wn))
        t5a = bool((wba - wbn) > cfg["t5a_min_delta"])
    else:
        wba = wbn = float("nan"); t5a = None

    out = dict(acc=float(acc), nir=float(nir), p_nir=float(p_nir), t1=t1,
               batch_only=float(bacc), t3_margin=float(acc - bacc), t3=t3,
               wb_acc=wba, wb_nir=wbn, t5a_delta=float(wba - wbn) if wb else None,
               t5a=t5a)
    if do_perm:
        r = _rng(seed + 500)
        null = np.array([cross_val_score(_clf(cfg), X, r.permutation(y),
                                         cv=_cv(cfg), n_jobs=1).mean()
                         for _ in range(n_perm)])
        out["perm_p95"] = float(np.percentile(null, 95))
        out["t2"] = bool(acc > out["perm_p95"])
    out["verdict"] = bool(out["t1"] and out["t3"] and (out["t5a"] is True))
    return out


def bh(p, q):
    p = np.asarray(p); o = np.argsort(p); m = p.size
    ok = p[o] <= q * np.arange(1, m + 1) / m
    out = np.zeros(m, bool)
    if ok.any():
        out[o[:np.max(np.where(ok)[0]) + 1]] = True
    return out


def cmd_audit(seeds, cfg=CONFIG):
    rows = []
    for s in seeds:
        for sc in ("A", "B"):
            X, y, b = simulate(sc, cfg["confound_default"], s, cfg=cfg)
            r = audit(X, y, b, cfg, do_perm=True, n_perm=25, seed=s)
            r.update(seed=s, scenario=sc)
            rows.append(r)
            print(f"seed={s} {sc}  acc={r['acc']:.3f} NIR={r['nir']:.3f} "
                  f"T1={'P' if r['t1'] else 'F'} T2={'P' if r.get('t2') else 'F'} "
                  f"T3margin={r['t3_margin']:+.3f}{'P' if r['t3'] else 'F'} "
                  f"T5delta={r['t5a_delta']:+.3f}{'P' if r['t5a'] else 'F'} "
                  f"-> {'PASS' if r['verdict'] else 'FAIL'}")
    return rows


def cmd_sweep(seeds, cfg=CONFIG):
    rows = []
    for c in [0.333, 0.50, 0.65, 0.80, 0.95]:
        for sc in ("A", "B"):
            accs, margins, deltas = [], [], []
            for s in seeds:
                X, y, b = simulate(sc, c, s, cfg=cfg)
                r = audit(X, y, b, cfg)
                accs.append(r["acc"]); margins.append(r["t3_margin"])
                deltas.append(r["t5a_delta"])
            rows.append(dict(confound=c, scenario=sc, n_seeds=len(seeds),
                             acc_mean=float(np.mean(accs)), acc_sd=float(np.std(accs)),
                             t3_margin_mean=float(np.mean(margins)),
                             t5a_delta_mean=float(np.mean(deltas)),
                             t5a_delta_sd=float(np.std(deltas))))
            print(f"conf={c:.3f} {sc}  acc={np.mean(accs):.3f}+-{np.std(accs):.3f} "
                  f"T3margin={np.mean(margins):+.3f}  "
                  f"T5delta={np.mean(deltas):+.3f}+-{np.std(deltas):.3f}")
    return rows


def cmd_floor(seeds, cfg=CONFIG):
    rows = []
    for eff in [0.10, 0.20, 0.40, 0.60, 0.80, 1.20]:
        sens, fdrs = [], []
        for s in seeds:
            X, y, b = simulate("A", cfg["confound_default"], s, effect=eff, cfg=cfg)
            pv = np.array([stats.kruskal(*[X[y == c, t] for c in range(cfg["n_cancers"])]).pvalue
                           for t in range(cfg["n_taxa"])])
            hit = bh(np.nan_to_num(pv, nan=1.0), cfg["fdr_q"])
            sig = np.arange(cfg["n_signal_taxa"])
            bat = np.arange(cfg["n_signal_taxa"], 2 * cfg["n_signal_taxa"])
            null = np.setdiff1d(np.arange(cfg["n_taxa"]), np.concatenate([sig, bat]))
            sens.append(hit[sig].mean()); fdrs.append(hit[null].mean())
        rows.append(dict(effect=eff, n_seeds=len(seeds),
                         sensitivity_mean=float(np.mean(sens)),
                         sensitivity_sd=float(np.std(sens)),
                         false_disc_mean=float(np.mean(fdrs))))
        print(f"effect={eff:.2f}  sens={np.mean(sens):.2f}+-{np.std(sens):.2f}  "
              f"falseDisc={np.mean(fdrs):.3f}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["audit", "sweep", "floor"])
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    print(f"# {CONFIG['version']}  cmd={a.cmd}  seeds={seeds}", file=sys.stderr)
    rows = {"audit": cmd_audit, "sweep": cmd_sweep, "floor": cmd_floor}[a.cmd](seeds)
    if a.out:
        json.dump({"config": CONFIG, "cmd": a.cmd, "seeds": seeds, "rows": rows},
                  open(a.out, "w"), indent=1)
