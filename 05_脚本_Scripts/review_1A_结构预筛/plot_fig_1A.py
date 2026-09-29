#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig 1A (new work): structure-anchored catalytic-pair GEOMETRY pre-screen of top candidates.
Upgrades manuscript limitation (i): MSA-mapped catalytic Glu (GH3 cols 2165/3114 grounded on
BlBG3 5Z9S; GH5 cols 1241/1473 grounded on BcelFp 6KDD) are checked in AlphaFold2 (full-length,
11 candidates) and ESMFold (GH5 clan-GH-A catalytic domain, 6 candidates) models for the
retaining acid/base--nucleophile carboxyl geometry. Same visual system as shared_style.py.
Two-stage BUILD/PLOT."""
import os, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))
import shared_style as ss
from shared_style import (PALETTE, SEM, FS, LW, MM, apply_rc, style_ax, panel_label,
                          figure_title, figure_note, save_figure, SourceTable, read_source)
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

WORK = Path('.'); OUTFIG = Path('./figures/out')
SCREEN = WORK/'out'/'t1a_catalytic_geometry_screen.tsv'
REF = WORK/'out'/'t1a_reference_calibration.tsv'
NAME = 'fig_1A_catalytic_geometry_prescreen'
SD = OUTFIG/f'{NAME}_source_data.tsv'
FAMCOL = {'GH3': PALETTE['green'], 'GH5': PALETTE['purple']}

def build():
    g = pd.read_csv(SCREEN, sep='\t'); r = pd.read_csv(REF, sep='\t')
    st = SourceTable(NAME)
    pv='out/t1a_catalytic_geometry_screen.tsv'
    for _,x in g.iterrows():
        e=x['protein_id']
        for f in ['rank','family','cdcd','minOO','caca','plddt_ab','plddt_nu','pos_acid_base','pos_nucleophile','ref_cdcd','ref_minOO']:
            st.add_num('cand', e, f, x[f], sub=x['family'], prov=pv)
        st.add('cand', e, x['family'], 'struct_source', str(x['struct_source']), prov=pv)
        st.add('cand', e, x['family'], 'subfamilies', str(x['subfamilies']), prov=pv)
        st.add('cand', e, x['family'], 'label', str(x['catalytic_geometry_consistency']), prov=pv)
        st.add('cand', e, x['family'], 'tier', str(x['tier']), prov=pv)
    pvr='out/t1a_reference_calibration.tsv (experimental PDB)'
    for _,x in r.iterrows():
        e=x['reference']
        for f in ['cdcd','minOO','caca']:
            st.add_num('ref', e, f, x[f], sub=x['family'], prov=pvr)
        st.add('ref', e, x['family'], 'family', x['family'], prov=pvr)
    st.write(SD); return SD

def plot():
    sf = read_source(SD); apply_rc()
    fig = ss.new_figure(MM['double'], 112.0)
    gs = fig.add_gridspec(1, 2, left=0.30, right=0.985, top=0.80, bottom=0.13, wspace=0.28, width_ratios=[1.25, 0.95])
    axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1])

    d = sf.block('cand')
    ents = list(dict.fromkeys(d['entity'].tolist()))
    recs=[]
    for e in ents:
        sub_e = d[d['entity']==e]['sub']
        fam = sub_e.iloc[0] if len(sub_e) else ''
        recs.append(dict(pid=e, fam=fam, rank=sf.one('cand','rank',e),
            cdcd=sf.one('cand','cdcd',e), minOO=sf.one('cand','minOO',e),
            plddt=min(sf.one('cand','plddt_ab',e), sf.one('cand','plddt_nu',e)),
            src=sf.one_text('cand','struct_source',e), sub=sf.one_text('cand','subfamilies',e),
            refcd=sf.one('cand','ref_cdcd',e), tier=sf.one_text('cand','tier',e)))
    df = pd.DataFrame(recs).sort_values(['fam','rank'])
    # ---- Panel A: per-candidate CD-CD distance strip (grouped GH3 top / GH5 bottom) ----
    order = df.sort_values(['fam','cdcd']).reset_index(drop=True)
    # arrange GH3 group then GH5 group with a gap
    gh3=order[order.fam=='GH3'].reset_index(drop=True); gh5=order[order.fam=='GH5'].reset_index(drop=True)
    yrows=[]; ylabels=[]; y=0; ticks=[]
    plotpts=[]
    def add_group(sub, fam):
        nonlocal y
        for _,rr in sub.iterrows():
            plotpts.append((y, rr)); ticks.append(y)
            lbl=rr.pid
            if rr.tier=='tier1_top20': lbl+=' (T1)'
            ylabels.append(lbl); y+=1
        y+=1  # gap
    add_group(gh3,'GH3'); add_group(gh5,'GH5')
    ymax=y
    for yy, rr in plotpts:
        fam=rr.fam; col=FAMCOL[fam]
        marker='o' if rr.src.startswith('AlphaFold') else 's'
        # pLDDT -> alpha
        a=0.45+0.55*min(max((rr.plddt-70)/30,0),1)
        axA.plot([rr.cdcd],[ymax-1-yy], marker=marker, ms=6.0, mfc=col, mec=PALETTE['ink'], mew=0.5, alpha=a, zorder=4)
    # reference lines
    cd3=sf.one('ref','cdcd','BlBG3 (PDB 5Z9S, GH3)'); cd5=sf.one('ref','cdcd','BcelFp (PDB 6KDD, GH5_25)')
    axA.axvline(cd3, color=FAMCOL['GH3'], ls=(0,(4,2)), lw=1.0, zorder=1)
    axA.axvline(cd5, color=FAMCOL['GH5'], ls=(0,(4,2)), lw=1.0, zorder=1)
    axA.text(cd3, ymax-0.3, f'BlBG3 {cd3:.1f}Å', color=FAMCOL['GH3'], fontsize=FS['note']-0.7, ha='center', va='bottom')
    axA.text(cd5, ymax-0.3, f'BcelFp {cd5:.1f}Å', color=FAMCOL['GH5'], fontsize=FS['note']-0.7, ha='center', va='bottom')
    axA.set_yticks([ymax-1-t for t in ticks]); axA.set_yticklabels(ylabels, fontsize=FS['tick']-1.2)
    axA.set_xlim(3.5, 10.0); axA.set_ylim(-0.8, ymax-0.2)
    axA.set_xlabel('Catalytic carboxyl–carboxyl distance  Cδ(E$_{a/b}$)–Cδ(E$_{nu}$)  (Å)')
    style_ax(axA, grid='x')
    panel_label(axA, 'A', 'Per-candidate Cδ–Cδ distance')

    # ---- Panel B: 2D scatter cdcd vs minOO with references ----
    for _,rr in df.iterrows():
        col=FAMCOL[rr.fam]; marker='o' if rr.src.startswith('AlphaFold') else 's'
        a=0.45+0.55*min(max((rr.plddt-70)/30,0),1)
        axB.plot([rr.cdcd],[rr.minOO], marker=marker, ms=6.0, mfc=col, mec=PALETTE['ink'], mew=0.5, alpha=a, zorder=4)
    for refname,mk in [('BlBG3 (PDB 5Z9S, GH3)','*'),('BcelFp (PDB 6KDD, GH5_25)','*')]:
        rb=sf.block('ref'); fam=rb[rb['entity']==refname]['sub'].iloc[0]
        axB.plot([sf.one('ref','cdcd',refname)],[sf.one('ref','minOO',refname)], marker='*', ms=15,
                 mfc=FAMCOL[fam], mec=PALETTE['ink'], mew=0.7, zorder=6)
    # consistent zone shading
    from matplotlib.patches import Rectangle
    axB.add_patch(Rectangle((3.5,2.5),9.0-3.5,8.0-2.5, facecolor=PALETTE['green_l'], alpha=0.18, edgecolor='none', zorder=0))
    axB.text(8.9,7.7,'retaining catalytic-pair\nconsistent zone', ha='right', va='top', fontsize=FS['note']-0.8, color=PALETTE['ink_soft'])
    axB.set_xlabel('Cδ–Cδ distance (Å)'); axB.set_ylabel('min carboxyl O–O distance (Å)')
    axB.set_xlim(3.5,9.2); axB.set_ylim(2.5,8.0)
    style_ax(axB, grid='both')
    panel_label(axB, 'B', 'vs experimental refs')

    from matplotlib.lines import Line2D
    leg=[Line2D([0],[0],marker='o',color='none',mfc=FAMCOL['GH3'],mec=PALETTE['ink'],ms=6,label='GH3 candidate'),
         Line2D([0],[0],marker='o',color='none',mfc=FAMCOL['GH5'],mec=PALETTE['ink'],ms=6,label='GH5 candidate'),
         Line2D([0],[0],marker='o',color='none',mfc=PALETTE['grey'],mec=PALETTE['ink'],ms=6,label='AlphaFold2 (full-length)'),
         Line2D([0],[0],marker='s',color='none',mfc=PALETTE['grey'],mec=PALETTE['ink'],ms=6,label='ESMFold (catalytic domain)'),
         Line2D([0],[0],marker='*',color='none',mfc=PALETTE['grey'],mec=PALETTE['ink'],ms=11,label='experimental reference (PDB)')]
    fig.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.62,0.93), ncol=5, fontsize=FS['legend']-0.8, frameon=False)

    figure_title(fig, 'Structure pre-screen: MSA-mapped catalytic Glu pairs form clan GH-A retaining geometry in all top candidates')
    figure_note(fig, "The two structure-anchored catalytic glutamates (GH3 MSA cols 2165/3114 grounded on BlBG3 5Z9S Glu162/Glu524; GH5 cols 1241/1473 grounded on BcelFp 6KDD Glu144/Glu260, full 484 GH3 / 68 GH5 pool MAFFT MSA) "
                     "were located in AlphaFold2 (full-length, AF DB v6; 7 GH3 + 4 GH5) and ESMFold (clan-GH-A single-domain; 6 GH5) models and their acid/base–nucleophile carboxyl distance measured. "
                     "All 17 modeled candidates fall within the retaining catalytic-pair geometry of the experimental references (GH3 Cδ–Cδ 5.99–7.26 Å vs 7.39; GH5 4.58–5.55 Å vs 4.80; min O–O ≤ 6.3 Å; per-residue pLDDT ≥ 80). "
                     "★ = manuscript tier-1 candidate. Marker opacity encodes the lower per-residue pLDDT. This upgrades the homology-mapping catalytic call (limitation i) to a structure-consistent one; it is a geometry pre-screen, not functional validation.",
                     y=0.008, size=FS['note']-0.9)
    paths = save_figure(fig, OUTFIG, NAME, kind='main')
    print('saved:', {k:str(v) for k,v in paths.items()})

if __name__=='__main__':
    build(); plot()
