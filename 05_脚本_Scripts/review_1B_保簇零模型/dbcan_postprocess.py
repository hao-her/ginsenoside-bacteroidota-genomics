#!/usr/bin/env python3
import os as _os9
_W9=_os9.environ.get('R9_WORK','work'); _I9=_os9.environ.get('R9_INPUT','input')
"""Post-processing: replicate dbcan process_results/filter_overlaps at a given threshold.
Input: raw scan TSV (from scan_proteome.py, domE=scan-level). Output: per-family unique-protein counts.
"""
import pandas as pd, numpy as np, sys

def filter_overlaps(df, thr=0.5):
    df=df.sort_values(['target_name','target_from','target_to'])
    kept=[]
    for _,g in df.groupby('target_name',sort=False):
        g=g.reset_index(drop=False).rename(columns={'index':'orig_idx'})
        if 'orig_idx' not in g.columns: g=g.reset_index(drop=True); g['orig_idx']=g.index
        keep=[]
        for i in range(len(g)):
            if not keep: keep.append(i); continue
            last=g.iloc[keep[-1]]; cur=g.iloc[i]
            overlap=min(last.target_to,cur.target_to)-max(last.target_from,cur.target_from)
            if overlap>0:
                len_last=max(1,last.target_to-last.target_from); len_cur=max(1,cur.target_to-cur.target_from)
                if overlap/len_last>thr or overlap/len_cur>thr:
                    if float(last.i_evalue)>float(cur.i_evalue): keep[-1]=i
                else: keep.append(i)
            else: keep.append(i)
        kept.extend(g.orig_idx.values[k] for k in keep)
    return df.loc[kept]

def counts_at(raw, thr):
    """raw: DataFrame with columns family,target_name,i_evalue,coverage... ; thr: i-evalue cutoff"""
    sub=raw[(raw.i_evalue<thr)&(raw.coverage>0.35)]
    if len(sub)==0: return {}
    f=sub.pipe(filter_overlaps)
    return f.groupby('family')['target_name'].nunique().to_dict()

if __name__=='__main__':
    raw=pd.read_csv(sys.argv[1],sep='\t')
    for t in [1e-15,1e-10,1e-5,1e-3]:
        c=counts_at(raw,t)
        print(f"thr={t:g}", {k:v for k,v in sorted(c.items()) if k.startswith('GH') and k in ['GH1','GH2','GH3','GH27','GH36','GH51','GH54','GH78','GH79','GH106']})
