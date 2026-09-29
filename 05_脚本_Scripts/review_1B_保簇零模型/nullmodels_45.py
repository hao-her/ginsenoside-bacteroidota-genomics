#!/usr/bin/env python3
"""Task 1B: cluster-preserving circular-shift null model (model C) on the UNBIASED
45-genome expanded set, alongside a model-P reproduction.

Windows are the delivered SusC-SusD all-pairs windows (Data/pul/t3_pairs_45.tsv);
genome-wide GH indicators are rebuilt from a fresh dbCAN.hmm V5-2 scan
(md5 8b277c30..., Z=875, i-evalue<1e-15, cov>0.35, cross-family overlap filter)
of the exact RefSeq proteomes, exactly matching the frozen pipeline.

Model P  (pair permutation; seed 20260921): reproduces Data/pul/t3_null45_aggregate.tsv.
Model C  (GH circular shift; seed 20260920, same family as r9_t2_condnull model C):
  windows FIXED; per replicon all GH ordinal positions shifted by a common uniform
  offset s in [0,L) (a bijection preserving every GH-GH spacing => GH clustering intact),
  scored via doubled-array prefix sums. NPERM=2000.
Empirical p = (1+#{null>=obs})/(1+NPERM). Pre-registered rule: ER>2 and p<0.01.
"""
import pandas as pd, numpy as np, sys, json, os
sys.path.insert(0,'.')
from dbcan_postprocess import filter_overlaps

WORK='.'; OUT=WORK+'/out'; os.makedirs(OUT,exist_ok=True)
GH10=sorted({f'GH{x}' for x in [1,2,3,27,36,51,54,78,79,106]})
CORE={'GH1','GH3'}; SIDE={'GH27','GH36','GH51','GH78','GH79','GH106'}
METRICS=['core','GH3','GH2','side','any']; NM=len(METRICS); NPERM=2000

P=pd.read_csv('原始数据包/Data/pul/t3_pairs_45.tsv',sep='\t')
r7=pd.read_csv('原始数据包/Data/matrices/official_strain_cazyme_matrix_round7.tsv',sep='\t',comment='#')
layer=dict(zip(r7.accession,r7.layer))
P['layer']=P.accession.map(layer)
accs=sorted(P.accession.unique())

# ---- build genome-wide GH indicators from fresh scan ----
ST={}
for acc in sorted(set(P.accession)|set(r7.accession)):
    genes=pd.read_csv(f'{WORK}/scan/gff/{acc}.genes.tsv',sep='\t')
    raw=pd.read_csv(f'{WORK}/scan/raw/{acc}.raw.tsv',sep='\t')
    f=filter_overlaps(raw[(raw.i_evalue<1e-15)&(raw.coverage>0.35)].copy())
    f=f[f.family.isin(GH10)]
    fam_of=f.groupby('target_name')['family'].apply(set).to_dict()
    pid2ord=dict(zip(genes.protein_id,genes.gene_ordinal))
    n=len(genes)
    ind=np.zeros((n,NM),dtype=bool)
    for pid,fams in fam_of.items():
        o=pid2ord.get(pid)
        if o is None: continue
        ind[o,0]=bool(fams&CORE); ind[o,1]='GH3' in fams; ind[o,2]='GH2' in fams
        ind[o,3]=bool(fams&SIDE); ind[o,4]=True
    ps=np.zeros((n+1,NM),dtype=np.int64); ps[1:]=np.cumsum(ind,axis=0)
    rr=genes.groupby('replicon')['gene_ordinal'].agg(['min','max']).to_dict('index')
    # doubled prefix per replicon for model C
    dbl={}
    for rep,r in rr.items():
        L=r['max']-r['min']+1
        sub=ind[r['min']:r['max']+1,:]
        d2=np.concatenate([sub,sub],axis=0)
        ps2=np.zeros((2*L+1,NM),dtype=np.int64); ps2[1:]=np.cumsum(d2,axis=0)
        dbl[rep]=(L,ps2)
    ST[acc]=dict(ind=ind,ps=ps,rr=rr,dbl=dbl)

# ---- FIDELITY GATE: recompute window GH content, compare to delivered ----
def win_families(acc,rep,w0,w1):
    ind=ST[acc]['ind']
    sub=ind[w0:w1+1,:]
    fams=[]
    # map columns to family strings via GH3/GH2/side detail requires per-gene family; recompute from raw
    return sub
# recompute has_core/has_gh/GH3/GH2 per window using prefix sums
def recompute_row(row):
    acc=row.accession; ps=ST[acc]['ps']; w0=int(row.w_start); w1=int(row.w_end)
    v=ps[w1+1]-ps[w0]
    return pd.Series(dict(rc_core=v[0]>0, rc_GH3=v[1]>0, rc_GH2=v[2]>0, rc_side=v[3]>0, rc_any=v[4]>0))
RC=P.apply(recompute_row,axis=1)
gate=dict(
  core=int((RC.rc_core==P.has_core).sum()),
  side=int((RC.rc_side==P.has_side).sum()),
  any=int((RC.rc_any==P.has_gh).sum()),
  GH3=int((RC.rc_GH3==P.gh_families.fillna('').str.split(';').apply(lambda l:'GH3' in l)).sum()),
  GH2=int((RC.rc_GH2==P.gh_families.fillna('').str.split(';').apply(lambda l:'GH2' in l)).sum()),
  n=len(P))
print("FIDELITY GATE (window recompute vs delivered):",gate,flush=True)
json.dump(gate,open(OUT+'/fidelity_gate_45.json','w'),indent=1)

# observed
obs={}
obs['any']=int(P.has_gh.sum()); obs['core']=int(P.has_core.sum()); obs['side']=int(P.has_side.sum())
obs['GH3']=int(P.gh_families.fillna('').str.split(';').apply(lambda l:'GH3' in l).sum())
obs['GH2']=int(P.gh_families.fillna('').str.split(';').apply(lambda l:'GH2' in l).sum())
print("observed (k15):",obs,"n_pairs:",len(P),flush=True)

# group windows
GRP={}
for (acc,rep),g in P.groupby(['accession','replicon']):
    lo0=np.minimum(g.o_susc.values,g.o_susd.values).astype(np.int64)
    d=(np.maximum(g.o_susc.values,g.o_susd.values)-lo0).astype(np.int64)
    GRP[(acc,rep)]=dict(rmin=ST[acc]['rr'][rep]['min'],rmax=ST[acc]['rr'][rep]['max'],
                        d=d,w0=g.w_start.values.astype(np.int64),w1=g.w_end.values.astype(np.int64),
                        wwidth=(g.w_end.values-g.w_start.values+1).astype(np.int64))

# ---- Model P (reproduce; seed 20260921) ----
rng=np.random.default_rng(20260921)
nullP={m:np.zeros(NPERM,dtype=np.int64) for m in METRICS}
nullP_layer={(l,m):np.zeros(NPERM,dtype=np.int64) for l in P.layer.unique() for m in ['core','any']}
for perm in range(NPERM):
    for (acc,rep),G in GRP.items():
        rmin,rmax,d=G['rmin'],G['rmax'],G['d']
        M=(rmax-d-rmin+1).astype(np.float64)
        lo=rmin+np.floor(rng.random(len(d))*M).astype(np.int64)
        w0=np.maximum(lo-10,rmin); w1=np.minimum(lo+d+10,rmax)
        ps=ST[acc]['ps']; has=(ps[w1+1]-ps[w0])>0
        lay=layer[acc]
        for mi,m in enumerate(METRICS): nullP[m][perm]+=has[:,mi].sum()
        nullP_layer[(lay,'core')][perm]+=has[:,0].sum(); nullP_layer[(lay,'any')][perm]+=has[:,4].sum()
    if (perm+1)%500==0: print("model P perm",perm+1,flush=True)

# ---- Model C (NEW; seed 20260920) ----
rng=np.random.default_rng(20260920)
nullC={m:np.zeros(NPERM,dtype=np.int64) for m in METRICS}
nullC_layer={(l,m):np.zeros(NPERM,dtype=np.int64) for l in P.layer.unique() for m in ['core','any']}
for perm in range(NPERM):
    for (acc,rep),G in GRP.items():
        L,ps2=ST[acc]['dbl'][rep]
        s=int(rng.integers(0,L))
        x0=G['w0']-G['rmin']; x0s=(x0-s)%L; w=G['wwidth']
        has=(ps2[x0s+w]-ps2[x0s])>0
        lay=layer[acc]
        for mi,m in enumerate(METRICS): nullC[m][perm]+=has[:,mi].sum()
        nullC_layer[(lay,'core')][perm]+=has[:,0].sum(); nullC_layer[(lay,'any')][perm]+=has[:,4].sum()
    if (perm+1)%500==0: print("model C perm",perm+1,flush=True)

def emp_p(o,nl): return (1+int((nl>=o).sum()))/(1+len(nl))
rows=[]
for model,null,seed in [('P_pairshuffle',nullP,20260921),('C_ghshift',nullC,20260920)]:
    for m in METRICS:
        nl=null[m]; o=obs[m]
        rows.append(dict(model=model,scope='all45',kou='k15',metric=m,observed=o,n_pairs=len(P),
                         null_mean=nl.mean(),null_sd=nl.std(),
                         null_p02_5=np.percentile(nl,2.5),null_p97_5=np.percentile(nl,97.5),
                         enrichment_ratio=o/nl.mean() if nl.mean()>0 else np.nan,
                         empirical_p_one_sided_ge=emp_p(o,nl),nperm=NPERM,seed=seed))
A=pd.DataFrame(rows)
A.to_csv(OUT+'/t3_null45_modelPC_aggregate.tsv',sep='\t',index=False)
pd.set_option('display.width',220); print(A.to_string())

# per-layer
lrows=[]
for model,nl_layer,seed in [('P_pairshuffle',nullP_layer,20260921),('C_ghshift',nullC_layer,20260920)]:
    for (l,m),nl in sorted(nl_layer.items()):
        sub=P[P.layer==l]; o=int(sub.has_core.sum()) if m=='core' else int(sub.has_gh.sum())
        lrows.append(dict(model=model,scope=l,kou='k15',metric=m,observed=o,n_pairs=len(sub),
                          null_mean=nl.mean(),null_sd=nl.std(),
                          enrichment_ratio=o/nl.mean() if nl.mean()>0 else np.nan,
                          empirical_p_one_sided_ge=emp_p(o,nl),nperm=NPERM,seed=seed))
L=pd.DataFrame(lrows); L.to_csv(OUT+'/t3_null45_modelPC_perlayer.tsv',sep='\t',index=False)
print(); print(L.to_string())
np.savez_compressed(OUT+'/t3_null45_modelC_distributions.npz',
    **{f"P|{m}":nullP[m] for m in METRICS},**{f"C|{m}":nullC[m] for m in METRICS},
    **{f"P|{l}|{m}":nullP_layer[(l,m)] for (l,m) in nullP_layer},
    **{f"C|{l}|{m}":nullC_layer[(l,m)] for (l,m) in nullC_layer})
print("\nsaved t3_null45_modelPC_aggregate.tsv, perlayer, distributions.npz")
