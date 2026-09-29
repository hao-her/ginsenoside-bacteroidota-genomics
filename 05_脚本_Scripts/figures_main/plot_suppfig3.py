#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_suppfig3.py — Suppl Fig 3 45 株无偏口径 PUL 逐株 / 逐层
=============================================================
【任务书要求】逐株配对占比 + 按目分层（L1/L2/L3）的富集对比；
标注 2 株 0 配对；无偏口径，不与条件化口径混比。

【设计思路】
Panel A：逐株点图。x = 45 株，按层（L1 Bacteroidales / L2 Flavobacteriales /
    L3 Cytophagales / L3 Sphingobacteriales）分组并铺层色带；
    y = 含 ≥1 GH 的 SusC–SusD 配对占比（大点），含 core 的占比（小点）；
    层均值以菱形标出。2 株 0 配对（无 SusC–SusD 配对）放入右侧独立 "no pairs" 槽，
    以空心灰点呈现 —— 其占比是"未定义"而非 0，绝不画成 0（红线 1）。
Panel B：逐层富集森林图（无偏全配对零模型，模型 P）。观测 ER 为点；零分布用高斯色彩晕
    （null_mean/null_sd）；有交付百分位者画 95% 置换区间须；ER=1 / ER=2 参考线；
    显著（ER>2 且 p<0.01）实心，否则空心灰。
注记：无偏口径 = 全 SusC–SusD 配对（不做 GH 预筛）；与 Fig 3A 的条件化口径不可混比（R7）。

运行：python plot_suppfig3.py / --plot-only
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import shared_style as S

FIG_NAME = "suppfig3_pul_unbiased_45strains"
OUT_DIR = S.SUPP_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 132.0
KIND = "main"

LAYERS = ["L1_Bacteroidales", "L2_Flavobacteriales", "L3_Cytophagales", "L3_Sphingobacteriales"]
LAYER_COLOR = {"L1_Bacteroidales": S.ORDER["Bacteroidales"],
               "L2_Flavobacteriales": S.ORDER["Flavobacteriales"],
               "L3_Cytophagales": S.ORDER["Cytophagales"],
               "L3_Sphingobacteriales": S.ORDER["Sphingobacteriales"]}


# ======================================================================
# BUILD
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Suppl Fig 3 — unbiased (all-pairs) PUL enrichment across 45 strains")
    P = {k: S.rel(S.SRC[k]) for k in ("strain45", "unbiased45", "unbiased45_layer", "pairs45")}

    ss = S.read_upstream("strain45")
    for _, r in ss.iterrows():
        st.add("A_strain", str(r["accession"]), str(r["layer"]), "n_pairs", str(r["n_pairs"]),
               text=f'{r["order"]}|{r["genus"]}', prov=P["strain45"])
        st.add("A_strain", str(r["accession"]), str(r["layer"]), "pairs_with_gh", str(r["pairs_with_gh"]), prov=P["strain45"])
        st.add("A_strain", str(r["accession"]), str(r["layer"]), "pairs_with_core", str(r["pairs_with_core"]), prov=P["strain45"])

    ag = S.read_upstream("unbiased45")
    for _, r in ag.iterrows():
        for f in ("observed", "n_pairs", "null_mean", "null_sd", "null_p02_5", "null_p97_5",
                  "enrichment_ratio", "empirical_p_one_sided_ge"):
            st.add("B_layer", f"all45_{r['metric']}", "all45", f, str(r[f]), prov=P["unbiased45"])
    ly = S.read_upstream("unbiased45_layer")
    for _, r in ly.iterrows():
        for f in ("observed", "n_pairs", "null_mean", "null_sd", "enrichment_ratio", "empirical_p_one_sided_ge"):
            st.add("B_layer", f"{r['scope']}_{r['metric']}", str(r["scope"]), f, str(r[f]), prov=P["unbiased45_layer"])
    st.add("meta", "er_ref", "no_enrichment", "value", "1.0", prov="definition")
    st.add("meta", "er_ref", "decision", "value", "2.0", prov="round8/9 criteria C3")
    for i, l in enumerate(LAYERS):
        st.add("meta", "layer_order", str(i), "name", text=l, prov="round9 t3 layer definition")

    st.write(SOURCE_DATA)
    print(f"[supp3] BUILD -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    A = sd.block("A_strain")
    B = sd.block("B_layer")
    er1 = sd.one("meta", "value", entity="er_ref", sub="no_enrichment")
    er2 = sd.one("meta", "value", entity="er_ref", sub="decision")
    layers = [t for t in (sd.one_text("meta", "name", entity="layer_order", sub=str(i))
                          for i in range(len(LAYERS))) if t]

    def aval(acc, layer, field):
        r = A[(A["entity"] == acc) & (A["sub"] == layer) & (A["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    def bval(key, scope, field):
        r = B[(B["entity"] == key) & (B["sub"] == scope) & (B["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    gsA = fig.add_gridspec(1, 1, left=0.055, right=0.560, top=0.905, bottom=0.335)
    axA = fig.add_subplot(gsA[0, 0])
    gsB = fig.add_gridspec(1, 1, left=0.700, right=0.985, top=0.905, bottom=0.335)
    axB = fig.add_subplot(gsB[0, 0])

    # ---------------- Panel A：逐株占比 ----------------
    x = 0.0
    xpos = {}
    zero_pair = []
    for L in layers:
        accs = A[(A["sub"] == L) & (A["field"] == "n_pairs")]["entity"].tolist()
        start = x
        for acc in accs:
            np_ = aval(acc, L, "n_pairs")
            gh = aval(acc, L, "pairs_with_gh")
            co = aval(acc, L, "pairs_with_core")
            if np_ == 0:
                zero_pair.append((acc, L))
                xpos[acc] = None
                continue
            xpos[acc] = x
            axA.scatter([x], [gh / np_], s=26, color=LAYER_COLOR[L], edgecolors="white",
                        linewidths=0.4, zorder=4)
            axA.scatter([x], [co / np_], s=9, color="white", edgecolors=LAYER_COLOR[L],
                        linewidths=0.8, zorder=5)
            x += 1.0
        # 层色带 + 层均值菱形
        ghfr, cofr, nn = [], [], 0
        for acc in accs:
            np_ = aval(acc, L, "n_pairs")
            if np_ and np_ > 0:
                ghfr.append(aval(acc, L, "pairs_with_gh") / np_)
                cofr.append(aval(acc, L, "pairs_with_core") / np_)
                nn += 1
        axA.axvspan(start - 0.5, x - 0.5, color=LAYER_COLOR[L], alpha=0.07, zorder=0)
        if nn:
            xc = (start + x - 1) / 2
            axA.scatter([xc], [np.mean(ghfr)], marker="D", s=52, color=LAYER_COLOR[L],
                        edgecolors=S.PALETTE["ink"], linewidths=0.7, zorder=6)
        axA.text((start + x - 1) / 2 + 0.4, -0.045, L.replace("_", " "), transform=axA.get_xaxis_transform(),
                 ha="right", va="top", fontsize=S.FS["note"], color=LAYER_COLOR[L],
                 fontweight="bold", rotation=32, rotation_mode="anchor")
        x += 1.2     # 层间空档
    # 0 配对槽
    gutter = x + 0.6
    for acc, L in zero_pair:
        axA.scatter([gutter], [0.0], s=30, facecolors="none", edgecolors=S.PALETTE["ink_faint"],
                    linewidths=0.9, zorder=5, marker="o")
    axA.text(gutter + 0.4, -0.045, "no pairs", transform=axA.get_xaxis_transform(), ha="right",
             va="top", fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], fontweight="bold",
             rotation=32, rotation_mode="anchor")
    axA.set_ylim(-0.06, 1.10)
    axA.set_xticks([])
    axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axA.set_ylabel("fraction of SusC–SusD pairs containing\n≥1 GH (large) / ≥1 core GH (small)",
                   fontsize=S.FS["axis"])
    S.style_ax(axA, grid="y", keep_spines=("bottom", "left"))
    S.panel_label(axA, "A", "Per-strain pair composition (45 genomes, unbiased all-pairs)")

    # ---------------- Panel B：逐层富集森林 ----------------
    rows = []
    for L in layers:
        for metric, col in (("core", LAYER_COLOR[L]), ("any", LAYER_COLOR[L])):
            rows.append((f"{L}_{metric}", L, metric, col))
    for metric in ("core", "GH3", "GH2", "side", "any"):
        rows.append((f"all45_{metric}", "all45", metric, S.SEM["unbiased"]))

    y = 0.0
    yticks, ylabels = [], []
    for key, scope, metric, col in rows:
        er = bval(key, scope, "enrichment_ratio")
        p = bval(key, scope, "empirical_p_one_sided_ge")
        nm = bval(key, scope, "null_mean"); nsd = bval(key, scope, "null_sd")
        p025, p975 = bval(key, scope, "null_p02_5"), bval(key, scope, "null_p97_5")
        sig = S.sig_mark(er, p)
        if np.isfinite(nsd) and nm > 0:
            S.gauss_halo(axB, y, 1.0, nsd / nm, col, span=3.0, height=0.60, alpha_peak=0.26,
                         x_clip=(0.0, 6.4))
        if np.isfinite(p025) and np.isfinite(p975) and nm > 0:
            axB.hlines(y, p025 / nm, p975 / nm, color=col, lw=0.9, alpha=0.85, zorder=3)
        elif np.isfinite(nsd) and nm > 0:
            axB.hlines(y, max(0, 1 - 1.96 * nsd / nm), 1 + 1.96 * nsd / nm, color=col,
                       lw=0.9, alpha=0.85, zorder=3)
        axB.scatter([er], [y], s=30, color=col if sig else "white", edgecolors=col,
                    linewidths=1.0, zorder=5, marker="o" if scope == "all45" else "D")
        yticks.append(y)
        ylabels.append(f"{scope.replace('L1_', '').replace('L2_', '').replace('L3_', '')} · {metric}")
        axB.text(6.5, y, S.pfmt(p), fontsize=S.FS["note"], va="center", ha="left",
                 color=col if sig else S.PALETTE["ink_faint"], fontweight="bold" if sig else "normal")
        y += 1.0
    axB.axvline(er1, color=S.PALETTE["ink_faint"], lw=0.8, zorder=1)
    axB.axvline(er2, color=S.PALETTE["ink_faint"], lw=0.8, ls=(0, (4, 2.5)), zorder=1)
    axB.set_yticks(yticks)
    axB.set_yticklabels(ylabels, fontsize=S.FS["tick"])
    axB.set_ylim(-0.8, y - 0.2)
    axB.invert_yaxis()
    axB.set_xlim(0, 6.4)
    axB.set_xticks([0, 1, 2, 3, 4, 5, 6])
    axB.set_xlabel("enrichment ratio (observed / null mean)", fontsize=S.FS["axis"])
    S.style_ax(axB, grid="x", keep_spines=("bottom", "left"))
    S.panel_label(axB, "B", "Per-layer enrichment (unbiased null, model P)")

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.055, right=0.985, top=0.185, bottom=0.010)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")
    np45 = int(bval("all45_core", "all45", "n_pairs"))
    axN.text(0, 1.0,
        f"Unbiased caliber: all {np45} SusC–SusD pairs of the 45 genomes are counted, with NO pre-filter on GH content — this is the primary PUL estimate and is NOT "
        "comparable with the conditional (GH-pre-filtered, 319-window) caliber of Fig 3A; the two answer different questions and are never pooled.",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", fontweight="bold")
    axN.text(0, 0.52,
        "Layer heterogeneity: Flavobacteriales shows the strongest core enrichment (ER 4.13), Bacteroidales other-genera and the two small L3 layers are not significant on core "
        "(low power), while 'any GH' is significant in L1/L2/L3-Sphingobacteriales.  Diamonds in panel A are layer means; open markers = core-GH fraction.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.10,
        f"{len(zero_pair)} of 45 strains yield no SusC–SusD pair at all; their fraction is undefined and they are shown in a separate gutter, not as zeros.  "
        "Sources: t3_strain_summary_45.tsv, t3_null45_aggregate.tsv, t3_null45_perlayer.tsv.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top")

    S.figure_title(fig, "Unbiased PUL–GH co-localisation across 45 Bacteroidota genomes and four orders",
                   x=0.055, y=0.988)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[supp3] PLOT  ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
