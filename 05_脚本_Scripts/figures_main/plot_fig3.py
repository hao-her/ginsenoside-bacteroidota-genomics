#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig3.py — Fig 3 PUL 基因组组织
====================================
【任务书要求】
(A) 三种零模型 + 无偏 vs 条件化的 enrichment ratio 森林图
(B) 属级分化散点/柱状图（Bacteroides vs Parabacteroides 的 GH51/GH106 比例 + Fisher p）
图注必须写明 "unbiased ER≈2.0 as primary estimate; conditional ER≈5.0 as upper bound"。

【设计思路】
Panel A：森林图，分两个**视觉隔离**的区块（红线 R7：口径不同不得混比）：
    上块 = 条件化口径（31 株 / 319 个"窗内必含 GH"的配对，round6 定义）：
           round8 均匀随机窗零模型、round9 模型 P（配对置换）、round9 模型 C（保簇循环移位）。
           → 这是"上界"。
    下块 = 无偏口径（全 SusC–SusD 配对，不做 GH 预筛）：31 株 1303 对 / 45 株 627 对。
           → 这是"主估计"。
    每行：观测 ER 为实心/空心点；零分布用**高斯色彩晕**呈现（null_mean/null_sd，
    任务书🟢鼓励的色彩晕替代传统误差棒）；95% 置换区间画细须；
    ER=1（无富集）与 ER=2（round8/9 criteria C3 预注册判定线）为参考线。
Panel B：蝴蝶图。Bacteroides（193 对）向左、Parabacteroides（126 对）向右，
    行 = 侧链家族；标注 OR 与 Fisher p（BH 校正）。GH78 作为"无分化"的阴性对照置灰。
    GH106 额外注明：分化是配对丰度差异而非有无差异（株级存在性 p 不显著）。

运行：python plot_fig3.py / --plot-only
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import shared_style as S

FIG_NAME = "fig3_pul_genomic_organization"
OUT_DIR = S.MAIN_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 148.0
KIND = "main"

GENUS_COL = {"Bacteroides": S.PALETTE["blue"], "Parabacteroides": S.PALETTE["purple"]}


# ======================================================================
# BUILD
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Fig 3 — PUL (SusC–SusD) genomic organisation of GH families")
    P = {k: S.rel(S.SRC[k]) for k in ("condnull", "uniformnull", "unbiased45",
                                      "unbiased31", "unbiased31_genus", "fisher_sidechain")}

    # --- A_forest：条件化口径（31 株 / 319 GH-filtered windows）---------
    un = S.read_upstream("uniformnull")
    r = un[(un["kou"] == "k15") & (un["metric"] == "core")].iloc[0]
    st.add("A_forest", "round8_uniform_window_null", "conditional", "er", str(r["enrichment_ratio"]),
           prov=P["uniformnull"])
    for f in ("observed", "null_mean", "null_sd", "null_p02_5", "null_p97_5", "empirical_p_one_sided_ge"):
        st.add("A_forest", "round8_uniform_window_null", "conditional", f, str(r[f]), prov=P["uniformnull"])
    st.add("A_forest", "round8_uniform_window_null", "conditional", "cohort", "",
           text="31 strains · 319 GH-filtered SusC–SusD windows (method-A definition)", prov=P["uniformnull"])

    cn = S.read_upstream("condnull")
    for model, label in (("P_pairshuffle", "round9_model_P_pair_shuffle"),
                         ("C_ghshift", "round9_model_C_cluster_preserving_shift")):
        rr = cn[(cn["model"] == model) & (cn["kou"] == "k15") & (cn["metric"] == "core")].iloc[0]
        st.add("A_forest", label, "conditional", "er", str(rr["enrichment_ratio"]), prov=P["condnull"])
        for f in ("observed", "null_mean", "null_sd", "null_p02_5", "null_p97_5", "empirical_p_one_sided_ge"):
            st.add("A_forest", label, "conditional", f, str(rr[f]), prov=P["condnull"])
        st.add("A_forest", label, "conditional", "cohort", "",
               text="31 strains · 319 GH-filtered SusC–SusD windows (method-A definition)", prov=P["condnull"])

    # --- A_forest：无偏口径（全配对）-----------------------------------
    u31 = S.read_upstream("unbiased31")
    rr = u31[(u31["kou"] == "k15") & (u31["metric"] == "core")].iloc[0]
    st.add("A_forest", "unbiased_31strains_allpairs", "unbiased", "er", str(rr["enrichment_ratio"]), prov=P["unbiased31"])
    for f in ("observed", "null_mean", "null_sd", "null_p02_5", "null_p97_5", "empirical_p_one_sided_ge"):
        st.add("A_forest", "unbiased_31strains_allpairs", "unbiased", f, str(rr[f]), prov=P["unbiased31"])
    st.add("A_forest", "unbiased_31strains_allpairs", "unbiased", "cohort", "",
           text="31 strains · 1303 all SusC–SusD pairs (no GH pre-filter)", prov=P["unbiased31"])

    u45 = S.read_upstream("unbiased45")
    for metric in ("core", "GH3", "GH2", "side", "any"):
        rr = u45[(u45["kou"] == "k15") & (u45["metric"] == metric)].iloc[0]
        st.add("A_forest", f"unbiased_45strains_allpairs_{metric}", "unbiased", "er",
               str(rr["enrichment_ratio"]), prov=P["unbiased45"])
        for f in ("observed", "null_mean", "null_sd", "null_p02_5", "null_p97_5", "empirical_p_one_sided_ge"):
            st.add("A_forest", f"unbiased_45strains_allpairs_{metric}", "unbiased", f, str(rr[f]), prov=P["unbiased45"])
        st.add("A_forest", f"unbiased_45strains_allpairs_{metric}", "unbiased", "cohort", "",
               text="45 strains · 627 all SusC–SusD pairs (no GH pre-filter)", prov=P["unbiased45"])

    ug = S.read_upstream("unbiased31_genus")
    for _, rr in ug.iterrows():
        if rr["metric"] != "core":
            continue
        st.add("A_forest", f"unbiased_31strains_{rr['scope']}", "unbiased", "er",
               str(rr["enrichment_ratio"]), prov=P["unbiased31_genus"])
        for f in ("observed", "null_mean", "null_sd", "empirical_p_one_sided_ge"):
            st.add("A_forest", f"unbiased_31strains_{rr['scope']}", "unbiased", f, str(rr[f]), prov=P["unbiased31_genus"])
        st.add("A_forest", f"unbiased_31strains_{rr['scope']}", "unbiased", "cohort", "",
               text=f"{rr['scope']} · all SusC–SusD pairs (n = {int(rr['n_pairs'])})", prov=P["unbiased31_genus"])

    # --- B_genus：配对级 + 株级存在性 Fisher ---------------------------
    fs = S.read_upstream("fisher_sidechain")
    for _, rr in fs.iterrows():
        lvl = "pair" if rr["level"] == "pair" else "strain_presence"
        e = f'{rr["family"]}'
        st.add("B_genus", e, lvl, "bact_n", str(rr["Bact_n"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "bact_d", str(rr["Bact_d"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "para_n", str(rr["Para_n"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "para_d", str(rr["Para_d"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "oddsratio", str(rr["oddsratio"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "p", str(rr["p"]), prov=P["fisher_sidechain"])
        st.add("B_genus", e, lvl, "p_BH", str(rr["p_BH"]), prov=P["fisher_sidechain"])

    # --- meta ----------------------------------------------------------
    st.add("meta", "er_ref", "no_enrichment", "value", "1.0", prov="definition")
    st.add("meta", "er_ref", "decision", "value", "2.0",
           text="round8/round9 criteria C3: ER>2 and p<0.01 => significant enrichment", prov="logs/criteria.md C3")
    for i, o in enumerate(["conditional", "unbiased"]):
        st.add("meta", "caliber_order", str(i), "name", text=o, prov="red line R7 separation")

    st.write(SOURCE_DATA)
    print(f"[fig3] BUILD  -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
PRETTY = {
    "round8_uniform_window_null": "round8 uniform random-window null",
    "round9_model_P_pair_shuffle": "round9 model P (SusC/SusD pair shuffle)",
    "round9_model_C_cluster_preserving_shift": "round9 model C (cluster-preserving GH shift)",
    "unbiased_31strains_allpairs": "31 strains, all pairs — core (GH1+GH3)",
    "unbiased_31strains_Bacteroides": "31 strains, all pairs — Bacteroides core",
    "unbiased_31strains_Parabacteroides": "31 strains, all pairs — Parabacteroides core",
    "unbiased_45strains_allpairs_core": "45 strains, all pairs — core",
    "unbiased_45strains_allpairs_GH3": "45 strains, all pairs — GH3",
    "unbiased_45strains_allpairs_GH2": "45 strains, all pairs — GH2",
    "unbiased_45strains_allpairs_side": "45 strains, all pairs — side-chain",
    "unbiased_45strains_allpairs_any": "45 strains, all pairs — any of 10 families",
}


def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    A = sd.block("A_forest").copy()
    B = sd.block("B_genus").copy()
    for df in (A, B):
        for c in ("value",):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    def aval(row_label, caliber, field):
        r = A[(A["entity"] == row_label) & (A["sub"] == caliber) & (A["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    er1 = sd.one("meta", "value", entity="er_ref", sub="no_enrichment")
    er2 = sd.one("meta", "value", entity="er_ref", sub="decision")

    cond_rows = [e for e in dict.fromkeys(A[(A["sub"] == "conditional")]["entity"])]
    unb_rows = [e for e in dict.fromkeys(A[(A["sub"] == "unbiased")]["entity"])]

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    gsA = fig.add_gridspec(1, 1, left=0.285, right=0.620, top=0.885, bottom=0.265)
    axA = fig.add_subplot(gsA[0, 0])
    gsB = fig.add_gridspec(1, 1, left=0.720, right=0.985, top=0.860, bottom=0.265)
    axB = fig.add_subplot(gsB[0, 0])

    # ---------------- Panel A：双口径森林图 ----------------
    rows = []
    y = 0.0
    # 条件化块
    block_top_cond = y
    for e in cond_rows:
        rows.append((y, e, "conditional")); y += 1.0
    y += 1.1                                   # 口径分隔空档
    block_top_unb = y
    for e in unb_rows:
        rows.append((y, e, "unbiased")); y += 1.0
    ymax = y

    # 区块底色（视觉隔离两口径）
    axA.axhspan(block_top_cond - 0.5, block_top_cond + len(cond_rows) - 0.5,
                color=S.SEM["conditional"], alpha=0.055, zorder=0)
    axA.axhspan(block_top_unb - 0.5, block_top_unb + len(unb_rows) - 0.5,
                color=S.SEM["unbiased"], alpha=0.055, zorder=0)

    axA.axvline(er1, color=S.PALETTE["ink_faint"], lw=0.8, ls="-", zorder=1)
    axA.axvline(er2, color=S.PALETTE["ink_faint"], lw=0.8, ls=(0, (4, 2.5)), zorder=1)
    axA.text(er1 - 0.08, ymax + 0.42, "ER = 1 (no enrichment)", fontsize=S.FS["note"],
             color=S.PALETTE["ink_faint"], ha="right", va="center")
    axA.text(er2 + 0.08, ymax + 0.42, "ER = 2 (pre-registered decision line)", fontsize=S.FS["note"],
             color=S.PALETTE["ink_faint"], ha="left", va="center")

    ylabels, ycolors = [], []
    for (yy, e, cal) in rows:
        col = S.SEM["conditional"] if cal == "conditional" else S.SEM["unbiased"]
        er = aval(e, cal, "er")
        nm = aval(e, cal, "null_mean"); nsd = aval(e, cal, "null_sd")
        p025, p975 = aval(e, cal, "null_p02_5"), aval(e, cal, "null_p97_5")
        p = aval(e, cal, "empirical_p_one_sided_ge")
        sig = S.sig_mark(er, p)

        # 高斯色彩晕：零分布在 ER 尺度上的密度脊（只用交付的 null_mean/null_sd）
        if np.isfinite(nsd) and nm > 0:
            S.gauss_halo(axA, yy, 1.0, nsd / nm, col, span=3.2, height=0.62,
                         alpha_peak=0.30, x_clip=(0.0, 8.2))
        # 95% 置换区间细须
        if np.isfinite(p025) and np.isfinite(p975) and nm > 0:
            axA.hlines(yy, p025 / nm, p975 / nm, color=col, lw=0.9, alpha=0.85, zorder=3)
            axA.vlines([p025 / nm, p975 / nm], yy - 0.13, yy + 0.13, color=col, lw=0.9, zorder=3)
        # 观测 ER
        axA.scatter([er], [yy], s=34, color=col if sig else "white",
                    edgecolors=col, linewidths=1.0, zorder=5,
                    marker="o" if cal == "conditional" else "D")
        ylabels.append(PRETTY.get(e, e))
        ycolors.append(col)
        # p 值直接标签
        axA.text(8.30, yy, S.pfmt(p), fontsize=S.FS["note"], va="center", ha="left",
                 color=col if sig else S.PALETTE["ink_faint"],
                 fontweight="bold" if sig else "normal")

    axA.set_yticks([r[0] for r in rows])
    axA.set_yticklabels(ylabels, fontsize=S.FS["tick"])
    for lbl, c in zip(axA.get_yticklabels(), ycolors):
        lbl.set_color(c)
    axA.set_ylim(-0.9, ymax + 0.75)
    axA.invert_yaxis()
    axA.set_xlim(0.0, 8.2)
    axA.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    axA.set_xlabel("enrichment ratio  (observed / null mean)", fontsize=S.FS["axis"])
    S.style_ax(axA, grid="x", keep_spines=("bottom", "left"))
    S.panel_label(axA, "A", "PUL enrichment — two calibers, kept apart")

    # 区块标题：横排放在两块之间的空档内（避免与行标签重叠）
    axA.text(0.06, block_top_cond - 0.62,
             "CONDITIONAL (ascertainment-biased) = UPPER BOUND",
             fontsize=S.FS["note"], color=S.SEM["conditional"], fontweight="bold",
             ha="left", va="center")
    axA.text(0.06, block_top_unb - 0.62,
             "UNBIASED (all pairs, no GH pre-filter) = PRIMARY ESTIMATE",
             fontsize=S.FS["note"], color=S.SEM["unbiased"], fontweight="bold",
             ha="left", va="center")

    # ---------------- Panel B：属级侧链分化（蝴蝶图）----------------
    def bval(fam, level, field):
        r = B[(B["entity"] == fam) & (B["sub"] == level) & (B["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    fams = ["GH51", "GH106", "GH78"]
    yb = np.arange(len(fams))[::-1] * 1.6
    for i, fam in enumerate(fams):
        bn, bd = bval(fam, "pair", "bact_n"), bval(fam, "pair", "bact_d")
        pn, pdn = bval(fam, "pair", "para_n"), bval(fam, "pair", "para_d")
        OR, pv, pbh = (bval(fam, "pair", k) for k in ("oddsratio", "p", "p_BH"))
        sig = pbh < 0.05
        fcol = S.FAMILY.get(fam, S.PALETTE["grey"])
        base = S.PALETTE["ink"] if sig else S.PALETTE["ink_faint"]

        # 左侧 Bacteroides，右侧 Parabacteroides
        axB.barh(yb[i], -bn / bd, height=0.62, color=GENUS_COL["Bacteroides"],
                 alpha=0.90 if sig else 0.35, edgecolor="white", linewidth=0.5, zorder=3)
        axB.barh(yb[i], pn / pdn, height=0.62, color=GENUS_COL["Parabacteroides"],
                 alpha=0.90 if sig else 0.35, edgecolor="white", linewidth=0.5, zorder=3)
        axB.text(-bn / bd - 0.012, yb[i], f"{int(bn)}/{int(bd)}", ha="right", va="center",
                 fontsize=S.FS["note"], color=GENUS_COL["Bacteroides"], fontweight="bold" if sig else "normal")
        axB.text(pn / pdn + 0.012, yb[i], f"{int(pn)}/{int(pdn)}", ha="left", va="center",
                 fontsize=S.FS["note"], color=GENUS_COL["Parabacteroides"], fontweight="bold" if sig else "normal")
        # 家族标签 + 角色色点
        axB.text(0.0, yb[i] + 0.46, fam, ha="center", va="bottom", fontsize=S.FS["axis"],
                 color=S.PALETTE["ink"], fontweight="bold")
        # OR / p 注记
        axB.text(0.0, yb[i] - 0.50,
                 f"OR = {OR:.2f}   {S.pfmt(pv)}" + (f"  (BH {S.pfmt(pbh)})" if sig else "  (n.s.)"),
                 ha="center", va="top", fontsize=S.FS["note"], color=base,
                 fontweight="bold" if sig else "normal")

    axB.axvline(0, color=S.PALETTE["ink_soft"], lw=0.8, zorder=2)
    axB.set_xlim(-0.30, 0.30)
    axB.set_ylim(-1.6, max(yb) + 1.2)
    axB.set_yticks([])
    axB.set_xticks([-0.25, -0.125, 0, 0.125, 0.25])
    axB.set_xticklabels(["0.25", "0.125", "0", "0.125", "0.25"])
    axB.set_xlabel("proportion of SusC–SusD pairs carrying the family\n(left: Bacteroides n=193 · right: Parabacteroides n=126)",
                   fontsize=S.FS["axis"])
    S.style_ax(axB, grid="x", keep_spines=("bottom",))
    S.panel_label(axB, "B", "Side-chain glycosidase divergence between the two genera")
    axB.text(-0.15, max(yb) + 1.05, "Bacteroides", ha="center", va="center",
             fontsize=S.FS["axis"], color=GENUS_COL["Bacteroides"], fontweight="bold")
    axB.text(+0.15, max(yb) + 1.05, "Parabacteroides", ha="center", va="center",
             fontsize=S.FS["axis"], color=GENUS_COL["Parabacteroides"], fontweight="bold")

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.055, right=0.985, top=0.185, bottom=0.012)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")

    u45c = aval("unbiased_45strains_allpairs_GH3", "unbiased", "er")
    cnd = aval("round9_model_C_cluster_preserving_shift", "conditional", "er")
    gh106_p = bval("GH106", "strain_presence", "p")
    gh51_p = bval("GH51", "strain_presence", "p")

    axN.text(0, 1.0,
        f"Unbiased ER ≈ 2.0 is the PRIMARY estimate; conditional ER ≈ 5.0 is an UPPER BOUND reflecting the ascertainment bias of method A "
        f"(its 319 windows are, by definition, GH-containing). The two calibers answer different questions and are never pooled (R7).",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", fontweight="bold")
    axN.text(0, 0.66,
        f"Gaussian colour-halos show the 2000-permutation null distribution in ER units (delivered null_mean / null_sd); thin whiskers = 95% permutation interval; "
        f"filled marker = significant under the pre-registered rule ER > 2 and p < 0.01. Observed: unbiased 45-strain GH3 ER = {u45c:.2f}; conditional model C core ER = {cnd:.2f}.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.32,
        f"Panel B: GH51 is enriched in Bacteroides and GH106 in Parabacteroides at the pair level (Fisher, BH-corrected). For GH106 the strain-presence test is not significant "
        f"(p = {gh106_p:.2f}), so the divergence is a pair-abundance difference, NOT a presence/absence difference; GH51 strain-presence p = {gh51_p:.1e} is only marginal. "
        "GH78 shows no genus-level divergence and is retained as a negative control. These are compositional / co-localisation evidence, not family-level functional claims.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")

    S.figure_title(fig, "PUL (SusC–SusD) organisation: GH families are non-randomly co-localised, and side-chain composition diverges between genera",
                   x=0.055, y=0.988, size=S.FS["fig_title"] - 0.5)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[fig3] PLOT   ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
