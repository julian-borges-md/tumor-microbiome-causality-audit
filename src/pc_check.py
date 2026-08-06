import numpy as np, pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
exec(open('wp0_tumor_vs_normal.py').read().split('if __name__')[0])
X, meta = load(); tn = taxnames()
name2id = {}
for k,v in tn.items():
    if isinstance(v,str): name2id.setdefault(v,k)
targets = {"COAD":["Fusobacterium nucleatum","Fusobacterium","Bacteroides fragilis","Escherichia coli","Streptococcus gallolyticus"],
           "STAD":["Helicobacter pylori"],
           "HNSC":["Fusobacterium nucleatum"]}
for proj, names in targets.items():
    m = meta["project"]==proj
    sub, subm = X.loc[m], meta.loc[m]
    tum = subm["sample_type"].str.contains("Tumor",na=False)
    nor = subm["sample_type"].str.contains("Normal",na=False)
    sel = tum|nor
    Xv = np.log1p(sub.loc[sel].values.astype(float)); y = tum[sel].astype(int).values
    cols = list(sub.columns)
    print(f"\n=== {proj}  tumor={y.sum()} normal={(1-y).sum()} ===")
    for nm in names:
        tid = name2id.get(nm)
        if tid is None: print(f"   {nm:32s} NOT IN TAXONOMY"); continue
        cand=[c for c in cols if str(c)==str(int(tid)) or str(c)==str(tid)]
        if not cand: print(f"   {nm:32s} taxid {int(tid)} absent from matrix"); continue
        j = cols.index(cand[0]); v = Xv[:,j]
        prev_t=(v[y==1]>0).mean(); prev_n=(v[y==0]>0).mean()
        if v.sum()==0: print(f"   {nm:32s} present in matrix but ZERO reads everywhere"); continue
        p=stats.mannwhitneyu(v[y==1],v[y==0],alternative='two-sided').pvalue
        print(f"   {nm:32s} prev tum={prev_t:.2f} norm={prev_n:.2f}  mean tum={v[y==1].mean():.3f} norm={v[y==0].mean():.3f}  p={p:.3f}")
# AUC version of tumor vs normal
print("\n=== tumor vs normal, AUC (imbalance-robust) ===")
for proj in ["COAD","STAD","HNSC","ESCA"]:
    m = meta["project"]==proj; sub, subm = X.loc[m], meta.loc[m]
    tum=subm["sample_type"].str.contains("Tumor",na=False); nor=subm["sample_type"].str.contains("Normal",na=False)
    sel=tum|nor; Xv=np.log1p(sub.loc[sel].values.astype(float)); Xv=Xv[:,Xv.sum(0)>0]
    y=tum[sel].astype(int).values
    if y.sum()<15 or (1-y).sum()<10: continue
    auc=cross_val_score(RandomForestClassifier(n_estimators=300,n_jobs=-1,random_state=0,min_samples_leaf=2),
        Xv,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),scoring='roc_auc',n_jobs=1).mean()
    print(f"   {proj}: AUC={auc:.3f}  (0.5 = no signal)   n_tum={y.sum()} n_norm={(1-y).sum()}")
