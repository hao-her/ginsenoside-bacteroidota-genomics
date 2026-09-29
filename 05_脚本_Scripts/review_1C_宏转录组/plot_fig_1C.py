#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig 1C (new work): public human-gut metatranscriptome expression evidence for the candidates.
473 passer proteins vs 25.2M predicted CDS from 3 human gut metatranscriptome samples
(EBI MGnify MGYS00005360, DIAMOND blastp). Same visual system as shared_style.py."""
import os, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))
import shared_style as ss
from shared_style import (PALETTE, SEM, FS, LW, MM, apply_rc, style_ax, panel_label,
                          figure_title, figure_note, save_figure, SourceTable, read_source)
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

WORK=Path('.'); OUTFIG=Path('./figures/out')
EXPR=WORK/'out'/'t1c_metatx_expression.tsv'
NAME='fig_1C_metatranscriptome_expression'
SD=OUTFIG/f'{NAME}_source_data.tsv'
CALL_COL={'expressed: near-exact transcript (>=95% id, >=80 aa)':PALETTE['blue'],
          'expressed: high-identity transcript (>=90% id, >=50 aa)':PALETTE['sky'],
          'homolog transcript detected (60-90% id)':PALETTE['orange'],
          'no transcript hit (>=60% id)':PALETTE['grey_l']}
CALL_SHORT={'expressed: near-exact transcript (>=95% id, >=80 aa)':'near-exact (≥95%)',
            'expressed: high-identity transcript (>=90% id, >=50 aa)':'high-identity (≥90%)',
            'homolog transcript detected (60-90% id)':'homolog (60–90%)',
            'no transcript hit (>=60% id)':'no hit'}
CALL_ORDER=list(CALL_COL.keys())
KEY=['WP_010801083.1','WP_258971837.1','WP_224264006.1','WP_195400649.1','WP_005860207.1',
     'WP_195474816.1','WP_005811975.1','WP_007216976.1','WP_258914935.1','WP_176553131.1']

def build():
    e=pd.read_csv(EXPR,sep='\t'); st=SourceTable(NAME); pv='out/t1c_metatx_expression.tsv'
    # Panel A: tier1 top20 + key candidates
    selA=e[(e.tier=='tier1_top20')|(e.protein_id.isin(KEY))].drop_duplicates('protein_id').sort_values('rank')
    for _,r in selA.iterrows():
        pid=r.protein_id
        for f in ['rank','best_pident','best_qcov','n_ge90_len50','n_ge95_len50','n_hits']:
            st.add_num('A',pid,f,r[f],sub=r['family'],prov=pv)
        st.add('A',pid,r['family'],'call','',text=str(r['metatx_expression_call']),prov=pv)
        st.add('A',pid,r['family'],'species','',text=str(r['species']),prov=pv)
    # Panel B: composition by family
    for fam in ['GH3','GH5']:
        sub=e[e.family==fam]
        for c in CALL_ORDER:
            st.add_num('B',fam,CALL_SHORT[c],int((sub.metatx_expression_call==c).sum()),sub=c,prov=pv)
    st.write(SD); return SD

def plot():
    sf=read_source(SD); apply_rc()
    fig=ss.new_figure(MM['double'],110.0)
    gs=fig.add_gridspec(1,2,left=0.20,right=0.985,top=0.80,bottom=0.14,wspace=0.30,width_ratios=[1.3,0.85])
    axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
    # Panel A
    dA=sf.block('A'); ents=list(dict.fromkeys(dA['entity'].tolist()))
    recs=[]
    for e in ents:
        sub=dA[dA.entity==e]; fam=sub['sub'].iloc[0]
        recs.append(dict(pid=e,fam=fam,rank=sf.one('A','rank',e),bp=sf.one('A','best_pident',e),
            n90=sf.one('A','n_ge90_len50',e),n95=sf.one('A','n_ge95_len50',e),
            call=sf.one_text('A','call',e)))
    df=pd.DataFrame(recs).sort_values('rank')
    y=np.arange(len(df))[::-1]
    for yy,(_,r) in zip(y,df.iterrows()):
        col=CALL_COL.get(r['call'],PALETTE['grey'])
        # size by n90 fragments
        sz=5+2.0*np.sqrt(max(r['n90'],0))
        axA.plot([r['bp']],[yy],marker='o',ms=min(sz,16),mfc=col,mec=PALETTE['ink'],mew=0.5,zorder=4)
    axA.axvline(95,color=PALETTE['blue'],ls=(0,(4,2)),lw=0.9,zorder=1)
    axA.axvline(90,color=PALETTE['sky'],ls=(0,(2,2)),lw=0.8,zorder=1)
    axA.text(95,len(df)-0.3,'95%',color=PALETTE['blue'],fontsize=FS['note']-0.6,ha='center',va='bottom')
    axA.text(90,len(df)-0.3,'90%',color=PALETTE['sky'],fontsize=FS['note']-0.6,ha='center',va='bottom')
    labels=[]
    for _,r in df.iterrows():
        lbl=f"{int(r['rank'])}. {r['pid']}"
        labels.append(lbl)
    axA.set_yticks(y); axA.set_yticklabels(labels,fontsize=FS['tick']-1.2)
    axA.set_xlim(60,102); axA.set_ylim(-0.8,len(df)-0.2)
    axA.set_xlabel('Best transcript-fragment identity in gut metatranscriptome (%)')
    style_ax(axA,grid='x')
    panel_label(axA,'A','Tier-1 & key candidates: expression evidence')
    # Panel B stacked
    fams=['GH3','GH5']; bottoms={f:0 for f in fams}
    xpos=np.arange(len(fams))
    for c in CALL_ORDER:
        vals=[sf.one('B',CALL_SHORT[c],fam,sub=c) for fam in fams]
        vals=[0 if np.isnan(v) else v for v in vals]
        axB.bar(xpos,vals,bottom=[bottoms[f] for f in fams],color=CALL_COL[c],edgecolor='white',lw=0.6,width=0.62,label=CALL_SHORT[c])
        for i,f in enumerate(fams): bottoms[f]+=vals[i]
    axB.set_xticks(xpos); axB.set_xticklabels([f"{f}\n(n={int(bottoms[f])})" for f in fams],fontsize=FS['tick'])
    axB.set_ylabel('Candidate passers'); style_ax(axB,grid='y')
    panel_label(axB,'B','All 473 passers by expression call')
    axB.legend(loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=2,fontsize=FS['legend']-1.0,frameon=False)

    figure_title(fig,'Public human-gut metatranscriptome supports in-vivo expression of candidate GH3/GH5 genes')
    figure_note(fig,"473 candidate passer proteins searched (DIAMOND blastp, --very-sensitive, e<1e-5, id≥60%) against 25.2M predicted CDS from 3 human faecal metatranscriptome samples "
                    "(EBI MGnify study MGYS00005360, BioProject PRJNA289586, pipeline 5.0). Metatranscriptome ORFs are short read-level fragments (~85 aa), so query coverage per hit is low (~10–12%); "
                    "identity over the aligned fragment and the number of independent high-identity fragments are the expression signal. 221/473 candidates have ≥90%-identity transcript evidence (87 near-exact ≥95%). "
                    "The high-identity Bacteroides candidates (WP_195474816 97.7%, WP_007216976 97.5%) show near-exact gut transcripts; the Parabacteroides tier-1 head shows only homolog-level hits in this cohort (lower genus abundance). "
                    "Marker size ∝ number of ≥90%-identity fragments. Absence of a hit is not evidence of non-expression (3 samples, one T1D cohort). This is a public-data expression proxy, not quantitative RNA-seq.",
                    y=0.008,size=FS['note']-0.9)
    paths=save_figure(fig,OUTFIG,NAME,kind='main'); print('saved:',{k:str(v) for k,v in paths.items()})

if __name__=='__main__':
    build(); plot()
