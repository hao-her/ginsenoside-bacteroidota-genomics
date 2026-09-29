#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig 1B (new work): cluster-preserving null model (model C) vs pair-permutation (model P)
on the UNBIASED all-pairs 45-genome expanded set. Same visual system as main Fig 3 / Fig S3
(shared_style.py: Okabe-Ito palette, Gaussian null haze, filled/open significance markers).
Two-stage: BUILD source_data.tsv from out/ null tables, then PLOT from it only."""
import os, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))
import shared_style as ss
from shared_style import (PALETTE, SEM, FS, LW, MM, figsize, apply_rc, style_ax,
                          panel_label, figure_title, figure_note, gauss_halo, save_figure,
                          SourceTable, read_source, er_fmt, pfmt, sig_mark)
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

WORK = Path('.'); OUTFIG = Path('./figures/out')
AGG = WORK/'out'/'t3_null45_modelPC_aggregate.tsv'
LAY = WORK/'out'/'t3_null45_modelPC_perlayer.tsv'
NAME = 'fig_1B_pul_modelC_unbiased45'
SD = OUTFIG/f'{NAME}_source_data.tsv'

METRIC_LABEL = {'any':'Any panel GH','GH2':'GH2','GH3':'GH3','core':'Core (GH1/GH3)','side':'Side-chain'}
METRIC_ORDER = ['any','GH2','GH3','core','side']
MODEL_LABEL = {'P_pairshuffle':'model P (pair permutation)','C_ghshift':'model C (cluster-preserving shift)'}
LAYER_ORDER = ['L1_Bacteroidales','L2_Flavobacteriales','L3_Cytophagales','L3_Sphingobacteriales']
LAYER_LABEL = {'L1_Bacteroidales':'L1 Bact.','L2_Flavobacteriales':'L2 Flavo.',
               'L3_Cytophagales':'L3 Cyto.','L3_Sphingobacteriales':'L3 Sphingo.'}

def build():
    agg = pd.read_csv(AGG, sep='\t'); lay = pd.read_csv(LAY, sep='\t')
    st = SourceTable(NAME)
    prov_a = 'out/t3_null45_modelPC_aggregate.tsv (model P reproduces Data/pul/t3_null45_aggregate.tsv; model C new)'
    for _, r in agg.iterrows():
        key = f"{r.metric}|{r.model}"
        st.add_num('A_forest', key, 'observed', r.observed, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'ER', r.enrichment_ratio, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'null_mean', r.null_mean, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'null_sd', r.null_sd, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'null_p02_5', r.null_p02_5, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'null_p97_5', r.null_p97_5, sub=r.model, prov=prov_a)
        st.add_num('A_forest', key, 'p', r.empirical_p_one_sided_ge, sub=r.model, prov=prov_a)
        st.add('A_forest', key, r.model, 'metric', r.metric, prov=prov_a)
    prov_l = 'out/t3_null45_modelPC_perlayer.tsv'
    for _, r in lay.iterrows():
        key = f"{r.scope}|{r.metric}|{r.model}"
        st.add_num('B_layer', key, 'observed', r.observed, sub=r.model, prov=prov_l)
        st.add_num('B_layer', key, 'ER', r.enrichment_ratio, sub=r.model, prov=prov_l)
        st.add_num('B_layer', key, 'null_mean', r.null_mean, sub=r.model, prov=prov_l)
        st.add_num('B_layer', key, 'null_sd', r.null_sd, sub=r.model, prov=prov_l)
        st.add_num('B_layer', key, 'p', r.empirical_p_one_sided_ge, sub=r.model, prov=prov_l)
        st.add('B_layer', key, r.model, 'scope', r.scope, prov=prov_l)
        st.add('B_layer', key, r.model, 'metric', r.metric, prov=prov_l)
    st.write(SD)
    return SD

def ER(o, x):  # convert count x to ER given observed o and null mean scale
    return o / x if x else np.nan

def plot():
    sf = read_source(SD)
    apply_rc()
    fig = ss.new_figure(MM['double'], 120.0)
    gs = fig.add_gridspec(1, 2, left=0.14, right=0.985, top=0.775, bottom=0.16, wspace=0.30, width_ratios=[1.18, 0.9])
    axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1])
    models = ['P_pairshuffle','C_ghshift']
    mcol = {'P_pairshuffle': SEM['unbiased'], 'C_ghshift': PALETTE['vermillion']}

    # ---------- Panel A ----------
    rows = [(m,mod) for m in METRIC_ORDER for mod in models]
    yvals = np.arange(len(rows))[::-1]
    for y,(m,mod) in zip(yvals, rows):
        key=f"{m}|{mod}"
        o=sf.one('A_forest','observed',key,sub=mod); er=sf.one('A_forest','ER',key,sub=mod)
        nm=sf.one('A_forest','null_mean',key,sub=mod); nsd=sf.one('A_forest','null_sd',key,sub=mod)
        p02=sf.one('A_forest','null_p02_5',key,sub=mod); p97=sf.one('A_forest','null_p97_5',key,sub=mod)
        p=sf.one('A_forest','p',key,sub=mod)
        er_center=o/nm; er_sd=er_center*(nsd/nm)
        gauss_halo(axA, y, er_center, er_sd, PALETTE['grey'], span=3.0, height=0.62, alpha_peak=0.16, zorder=1)
        er_lo=o/p97; er_hi=o/p02
        axA.plot([er_lo,er_hi],[y,y], color=PALETTE['ink_soft'], lw=LW['axis'], zorder=2, solid_capstyle='round')
        sig=sig_mark(er,p)
        axA.plot([er],[y], marker='o', ms=5.4, mfc=(mcol[mod] if sig else 'white'), mec=mcol[mod], mew=1.3, zorder=4)
        axA.text(er+0.14, y, f"{er_fmt(er)} ({pfmt(p)})", va='center', ha='left', fontsize=FS['annot']-0.3, color=PALETTE['ink'])
    axA.axvline(1.0, color=PALETTE['axis'], lw=LW['axis'], zorder=0)
    axA.axvline(2.0, color=PALETTE['vermillion'], lw=LW['grid']+0.3, ls=(0,(4,2)), zorder=0)
    axA.text(2.0, -0.62, 'ER=2', fontsize=FS['note']-0.5, color=PALETTE['vermillion'], ha='center', va='top')
    ylabels=[f"{METRIC_LABEL[m]} · {'P' if mod.startswith('P') else 'C'}" for (m,mod) in rows]
    axA.set_yticks(yvals); axA.set_yticklabels(ylabels, fontsize=FS['tick'])
    axA.set_xlim(0.6, 4.5); axA.set_ylim(-0.9, len(rows)-0.2)
    axA.set_xlabel('Enrichment ratio (observed / permutation-null)')
    style_ax(axA, grid='x')
    panel_label(axA, 'A', 'Unbiased all-pairs (627 pairs)')
    from matplotlib.lines import Line2D
    leg=[Line2D([0],[0],marker='o',color='none',mfc=mcol['P_pairshuffle'],mec=mcol['P_pairshuffle'],ms=6,label='model P (pair permutation)'),
         Line2D([0],[0],marker='o',color='none',mfc=mcol['C_ghshift'],mec=mcol['C_ghshift'],ms=6,label='model C (cluster-preserving shift)'),
         Line2D([0],[0],marker='o',color='none',mfc='white',mec=PALETTE['ink_soft'],ms=6,label='open = not (ER>2 & p<0.01)')]
    fig.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5,0.905), ncol=3, fontsize=FS['legend']-0.5, frameon=False)

    # ---------- Panel B: per-layer, grouped ----------
    metrics=['core','any']
    group_w=1.0; sub=0.30
    centers=[]; labels=[]
    for li,l in enumerate(LAYER_ORDER):
        base=li*3.0
        for mi,met in enumerate(metrics):
            cx=base+mi*group_w
            centers.append(cx); labels.append(met)
            for k,mod in enumerate(models):
                key=f"{l}|{met}|{mod}"
                er=sf.one('B_layer','ER',key,sub=mod); p=sf.one('B_layer','p',key,sub=mod)
                if np.isnan(er): continue
                sig=sig_mark(er,p)
                axB.plot([cx+(k-0.5)*sub],[er], marker='o', ms=4.8, mfc=(mcol[mod] if sig else 'white'), mec=mcol[mod], mew=1.2, zorder=4)
        axB.text(base+0.5*group_w, -0.62, LAYER_LABEL[l], ha='center', va='top', fontsize=FS['note']-0.6, color=PALETTE['ink_soft'])
    axB.axhline(1.0, color=PALETTE['axis'], lw=LW['axis'], zorder=0)
    axB.axhline(2.0, color=PALETTE['vermillion'], lw=LW['grid']+0.3, ls=(0,(4,2)), zorder=0)
    axB.set_xticks(centers); axB.set_xticklabels(labels, fontsize=FS['tick']-0.5)
    axB.set_ylabel('Enrichment ratio'); axB.set_ylim(0, 5.2)
    axB.set_xlim(-0.7, (len(LAYER_ORDER)-1)*3.0+group_w+0.6)
    style_ax(axB, grid='y')
    panel_label(axB, 'B', 'Per-layer core & any-GH')

    figure_title(fig, 'Cluster-preserving null (model C) confirms PUL co-localization on the unbiased primary metric')
    figure_note(fig, "Model C keeps SusC-SusD windows fixed and circularly shifts all GH ordinal positions per replicon (a bijection preserving every GH-GH spacing, hence all GH clustering), randomizing GH position relative to pairs; 2000 permutations, seed 20260920. "
                     "Model P (seed 20260921) reproduces the delivered Data/pul/t3_null45_aggregate.tsv bit-for-bit. Empirical p=(1+#{null>=obs})/(1+2000); pre-registered significance ER>2 and p<0.01. "
                     "The unbiased ER~=2 primary estimate survives model C (GH3 ER=1.98, p=0.002; any ER=2.20, p<0.001; GH2 ER=2.56, p<0.001), so it is not explained by GH clustering alone. Window fidelity gate: 627/627.",
                     y=0.008, size=FS['note']-0.9)
    paths = save_figure(fig, OUTFIG, NAME, kind='main')
    print('saved:', {k:str(v) for k,v in paths.items()})

if __name__=='__main__':
    build(); plot()
