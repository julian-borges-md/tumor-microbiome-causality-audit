import numpy as np, pandas as pd
from scipy import stats
exec(open('wp0_nesting.py').read().split('if __name__')[0])
X, meta = load(); tax = tax_table()

def sig_set(proj, paired=True):
    if paired:
        T,N,cols = paired_matrix(X, meta, proj)
        if T is None or T.shape[0]<5: return None,None,None
        nz=(T.sum(0)+N.sum(0))>0; Ts,Ns,ids=T[:,nz],N[:,nz],cols[nz]
        pv=np.ones(Ts.shape[1])
        for j in range(Ts.shape[1]):
            d=Ts[:,j]-Ns[:,j]
            if np.count_nonzero(d)<3: continue
            try: pv[j]=stats.wilcoxon(d,zero_method="wilcox").pvalue
            except Exception: pv[j]=1.0
        M=np.vstack([Ts,Ns])
    else:
        m=meta["project"]==proj; sub,subm=X.loc[m],meta.loc[m]
        tum=subm["sample_type"].str.contains("Tumor",na=False)
        nor=subm["sample_type"].str.contains("Normal",na=False)
        sel=tum|nor; Xv=np.log1p(sub.loc[sel].values.astype(float))
        kc=Xv.sum(0)>0; Xv,ids=Xv[:,kc],np.array(sub.columns)[kc]
        y=tum[sel].astype(int).values
        pv=np.array([stats.mannwhitneyu(Xv[y==1,j],Xv[y==0,j],alternative="two-sided").pvalue for j in range(Xv.shape[1])])
        M=Xv
    pv=np.nan_to_num(pv,nan=1.0); hit=bh(pv,0.05)
    return list(ids[hit]), M, list(ids)

print("=== DISCOVERY COUNT vs REDUNDANCY RULE ===")
print(f"{'analysis':<22} {'reported':>9} " + " ".join(f"{f'r>{t}':>8}" for t in [0.99,0.95,0.90,0.80]) + f"{'ancestor':>10}")
for label,proj,paired in [("STAD paired","STAD",True),("HNSC unpaired","HNSC",False)]:
    hids,M,allids = sig_set(proj,paired)
    if not hids: continue
    sub_idx=[allids.index(h) for h in hids]; Mh=M[:,sub_idx]
    row=[]
    for t in [0.99,0.95,0.90,0.80]:
        g,_=redundancy_graph(hids,tax,Mh,rthresh=t); row.append(len(g))
    # pure ancestor rule
    lin={i:(tax.loc[float(i),"lineage"] if float(i) in tax.index else str(i)) for i in hids}
    cl=[]
    for i in sorted(hids,key=lambda z:len(lin[z])):
        p=False
        for c in cl:
            if lin[i].startswith(lin[c[0]]) or lin[c[0]].startswith(lin[i]): c.append(i);p=True;break
        if not p: cl.append([i])
    print(f"{label:<22} {len(hids):>9} " + " ".join(f"{v:>8}" for v in row) + f"{len(cl):>10}")

print("\n=== FEATURE-SPACE REDUNDANCY (paired cohorts) ===")
print(f"{'cancer':>7} {'nominal':>8} " + " ".join(f"{f'r>{t}':>8}" for t in [0.99,0.95,0.90]))
for proj in ["COAD","STAD","HNSC","ESCA"]:
    T,N,cols=paired_matrix(X,meta,proj)
    if T is None: continue
    nz=(T.sum(0)+N.sum(0))>0; ids=list(cols[nz]); M=np.vstack([T[:,nz],N[:,nz]])
    if len(ids)>500: ids=ids[:500]; M=M[:,:500]
    row=[]
    for t in [0.99,0.95,0.90]:
        g,_=redundancy_graph(ids,tax,M,rthresh=t); row.append(len(g))
    print(f"{proj:>7} {len(ids):>8} " + " ".join(f"{v:>8}" for v in row))
