#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig1.py — Fig 1 门级 / 属级 GH 家族拷贝数分布
====================================================
【任务书要求】箱线图，横轴=门（或属），纵轴=拷贝数（log 轴），
分别展示 GH1 / GH2 / GH3（三联图）；不标显著性星号；
图注标注 "descriptive pattern + strong phylogenetic signal (D≈0.02, λ≈0.93–0.99)"。

【设计思路】
Panel A/B/C（上排三联）：门级。箱体几何**直接取自已交付的**
    phylum_family_distribution.tsv 的 q1/median/q3/min/max（零重算），
    叠加 official_strain_cazyme_matrix_round4.tsv 的株级原始点（jitter）。
    Bacteroidota 用醒目描边 + 底色带强调 —— 这是本文核心谱系。
Panel D（下排通栏）：属级 12 属哑铃图，三个家族各一个 marker，按门着色分组。
    这一面板存在的理由：正是"134 株折叠成 12 属中位数"这一步让门级
    Kruskal–Wallis 失去显著性；把它画出来读者才能理解为何不标星号。
    数据取自 t3a_genus_medians.tsv。

纵轴用 symlog（linthresh=1）而非纯 log：拷贝数含大量真 0（GH1 在
Bacteroidota 31/31 株为 0），纯 log 轴会把它们丢掉 —— 那是选择性呈现，
违反红线 1。symlog 在 [0,1] 段线性、其上对数，0 与 1 都被如实画出，
轴标签明确写出该事实。

运行：python plot_fig1.py            （BUILD + PLOT）
      python plot_fig1.py --plot-only（仅从 source_data.tsv 重绘）
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import shared_style as S

# ----------------------------------------------------------------------
# 路径与版面变量（顶部定义，不硬编码进函数）
# ----------------------------------------------------------------------
FIG_NAME = "fig1_phylum_genus_gh_distribution"
OUT_DIR = S.MAIN_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 150.0     # 双栏 183 mm
KIND = "main"                                    # 主图 → 600 dpi

FAMILIES = ["GH1", "GH2", "GH3"]
MARKER = {"GH1": "o", "GH2": "s", "GH3": "D"}
JITTER_SEED = 20260920       # 固定种子：jitter 可复现，且不改变任何数值


# ======================================================================
# 阶段 1 —— BUILD：从上游只读表装配 source_data.tsv
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Fig 1 — phylum/genus GH family copy-number distribution")

    m4 = S.read_upstream("matrix_r4")
    ph = S.read_upstream("phylum_dist")
    ge = S.read_upstream("genus_dist")
    gm = S.read_upstream("genus_medians")
    gkw = S.read_upstream("genus_kw")
    skw = S.read_upstream("strain_kw")
    dstat = S.read_upstream("phylo_D")
    lam = S.read_upstream("phylo_lambda")

    P_PH, P_M4, P_GM = S.rel(S.SRC["phylum_dist"]), S.rel(S.SRC["matrix_r4"]), S.rel(S.SRC["genus_medians"])

    # --- phylum_box：门级箱体几何（逐字搬运已交付汇总表，不重算）-------
    for _, r in ph.iterrows():
        if r["family"] not in FAMILIES:
            continue
        for f in ("n", "median", "q1", "q3", "iqr", "min", "max", "mean", "zero_frac"):
            st.add("phylum_box", r["phylum"], r["family"], f, str(r[f]), prov=P_PH)

    # --- phylum_points：株级原始值（jitter 叠加层）---------------------
    for _, r in m4.iterrows():
        for fam in FAMILIES:
            st.add("phylum_points", str(r["accession"]), fam, "value", str(r[fam]),
                   text=f'{r["genus"]}|{r["phylum"]}', prov=P_M4)

    # --- genus_median / genus_n：属级哑铃图 ----------------------------
    for _, r in gm.iterrows():
        for fam in FAMILIES:
            st.add("genus_median", str(r["genus"]), fam, "median", str(r[fam]),
                   text=str(r["phylum"]), prov=P_GM)
    for _, r in ge.iterrows():
        if r["family"] in FAMILIES:
            st.add("genus_n", str(r["genus"]), r["family"], "n", str(r["n"]),
                   prov=S.rel(S.SRC["genus_dist"]))

    # --- stats：图注所需的已定稿统计量（只搬运，不计算）----------------
    for _, r in gkw.iterrows():
        if r["family"] in FAMILIES:
            st.add("stats", "genus_median_phylum_KW", str(r["family"]), "p", str(r["p"]),
                   text=f'H={r["H"]}; k={r["k"]}; N={r["N"]}', prov=S.rel(S.SRC["genus_kw"]))
    for _, r in skw.iterrows():
        if r["family"] in FAMILIES and r["level"] == "phylum":
            st.add("stats", "strain_phylum_KW", str(r["family"]), "epsilon_squared",
                   str(r["epsilon_squared"]), text=f'H={r["H"]}; p={r["p_value"]}',
                   prov=S.rel(S.SRC["strain_kw"]))
    for _, r in dstat.iterrows():
        st.add("stats", "phylo_signal", str(r["trait"]), "D", str(r["D"]),
               text=f'p_bm={r["p_bm"]}; p_rand={r["p_rand"]}', prov=S.rel(S.SRC["phylo_D"]))
    for _, r in lam.iterrows():
        st.add("stats", "phylo_signal", str(r["family"]), "lambda_ML", str(r["lambda_ML"]),
               text=f'CI={r["lambda_CI_lo"]}-{r["lambda_CI_hi"]}; K={r["K"]}; K_perm_p={r["K_perm_p"]}',
               prov=S.rel(S.SRC["phylo_lambda"]))

    # --- meta：版面参数（非数据）--------------------------------------
    for i, p in enumerate(S.PHYLUM_ORDER):
        st.add("meta", entity="phylum_order", sub=str(i), field="name", text=p,
               prov="shared_style.PHYLUM_ORDER")
    for i, f in enumerate(FAMILIES):
        st.add("meta", entity="family_order", sub=str(i), field="name", text=f,
               prov="task spec (GH1/GH2/GH3)")
    st.add("meta", entity="y_axis", field="scale", value="1.0", text="symlog",
           prov="design decision: linthresh=1 keeps true zeros visible, not dropped")
    st.add("meta", entity="jitter_seed", field="value", value=str(JITTER_SEED),
           prov="reproducibility")

    st.write(SOURCE_DATA)
    print(f"[fig1] BUILD  -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# 阶段 2 —— PLOT：只从 source_data.tsv 读数
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)

    phyla = [t for t in (sd.one_text("meta", "name", entity="phylum_order", sub=str(i))
                         for i in range(len(S.PHYLUM_ORDER))) if t]
    fams = [t for t in (sd.one_text("meta", "name", entity="family_order", sub=str(i))
                        for i in range(len(FAMILIES))) if t]
    seed = int(sd.one("meta", "value", entity="jitter_seed"))
    yscale = sd.one_text("meta", "scale", entity="y_axis")
    ylinthresh = sd.one("meta", "scale", entity="y_axis")

    A = sd.block("phylum_box").copy()
    Pts = sd.block("phylum_points").copy()
    Bm = sd.block("genus_median").copy()
    Bn = sd.block("genus_n").copy()
    for _df in (A, Pts, Bm):
        _df["v"] = pd.to_numeric(_df["value"], errors="coerce")

    def boxval(ph, fam, field):
        r = A[(A["entity"] == ph) & (A["sub"] == fam) & (A["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    rng = np.random.default_rng(seed)
    hl = phyla.index("Bacteroidota") if "Bacteroidota" in phyla else 0

    # ---------------- 版面：三行 gridspec，第三行专放强制注记，保证不压图 ----
    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    # 两个独立 gridspec：绘图区 + 注记区（注记绝不压图）
    gsA = fig.add_gridspec(2, 3, left=0.052, right=0.995, top=0.905, bottom=0.235,
                           hspace=0.46, wspace=0.20, height_ratios=[1.0, 1.15])
    axA = [fig.add_subplot(gsA[0, j]) for j in range(3)]
    axB = fig.add_subplot(gsA[1, :])
    gsN = fig.add_gridspec(1, 1, left=0.052, right=0.995, top=0.155, bottom=0.012)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")

    # ============ Panel A/B/C：门级箱线 + 株级原始点 ============
    for j, fam in enumerate(fams):
        ax = axA[j]
        S.band(ax, hl - 0.5, hl + 0.5, S.SEM["not_detected"], alpha=0.075)

        sub = Pts[Pts["sub"] == fam]
        for i, ph in enumerate(phyla):
            col = S.PHYLUM.get(ph, S.PALETTE["grey"])
            q1, med, q3 = (boxval(ph, fam, k) for k in ("q1", "median", "q3"))
            lo, hi = boxval(ph, fam, "min"), boxval(ph, fam, "max")
            nn = int(boxval(ph, fam, "n"))
            emph = (i == hl)
            bw = 0.46 if emph else 0.40

            g = sub[sub["text"].str.split("|").str[1] == ph]["v"].to_numpy(dtype=float)
            if len(g):
                ax.scatter(i + rng.uniform(-0.20, 0.20, len(g)), g,
                           s=6.5 if emph else 5.0, color=col,
                           alpha=0.62 if emph else 0.48,
                           edgecolors="white", linewidths=0.22, zorder=3)

            ax.vlines(i, lo, hi, color=col, lw=0.85, alpha=0.85, zorder=4)
            ax.hlines([lo, hi], i - bw * 0.30, i + bw * 0.30, color=col,
                      lw=0.85, alpha=0.85, zorder=4)
            ax.add_patch(Rectangle((i - bw / 2, q1), bw, max(q3 - q1, 1e-9),
                                   facecolor=col, edgecolor=S.PALETTE["ink"] if emph else "white",
                                   linewidth=1.05 if emph else 0.7,
                                   alpha=0.92 if emph else 0.76, zorder=5))
            ax.hlines(med, i - bw / 2, i + bw / 2,
                      color=S.PALETTE["ink"] if emph else "white",
                      lw=1.5 if emph else 1.25, zorder=6, capstyle="butt")
            ax.text(i + bw / 2 + 0.06, med, f"{med:g}", fontsize=S.FS["note"],
                    color=S.PALETTE["ink"] if emph else S.PALETTE["ink_soft"],
                    va="center", ha="left", zorder=7,
                    fontweight="bold" if emph else "normal")
            # 样本量放在面板顶部（避免与旋转的门名打架）
            ax.text(i, 0.985, f"n={nn}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=S.FS["note"],
                    color=S.PALETTE["ink_faint"], zorder=8,
                    bbox=dict(boxstyle="round,pad=0.12", fc=S.PALETTE["paper"],
                              ec="none", alpha=0.85))

        ax.set_xticks(range(len(phyla)))
        ax.set_xticklabels(phyla, fontsize=S.FS["tick"], rotation=26, ha="right")
        for lbl, p in zip(ax.get_xticklabels(), phyla):
            lbl.set_color(S.PHYLUM.get(p, S.PALETTE["ink_soft"]))
            if p == "Bacteroidota":
                lbl.set_fontweight("bold")
        ax.set_xlim(-0.62, len(phyla) - 0.20)
        ax.set_yscale(yscale, linthresh=ylinthresh)
        ax.set_yticks([0, 1, 2, 5, 10, 20, 50])
        ax.set_ylim(-1.4, 62)
        S.style_ax(ax, grid="y")
        S.panel_label(ax, "ABC"[j], f"{fam}  \u00b7  {S.FAMILY_ROLE[fam]} family")
        if j == 0:
            ax.set_ylabel("GH copy number per genome\n(symlog: linear 0\u20131, log above)",
                          fontsize=S.FS["axis"])
        ax.set_xlabel("Phylum", fontsize=S.FS["axis"], labelpad=1.5)

    # ============ Panel D：属级哑铃（12 属中位数）============
    ax = axB
    genera = []
    for ph in phyla:
        gg = Bm[Bm["text"] == ph]["entity"].drop_duplicates().tolist()
        gg.sort(key=lambda g: -float(Bm[(Bm["entity"] == g) & (Bm["sub"] == "GH2")]["v"].iloc[0]))
        genera += [(g, ph) for g in gg]

    ypos = np.arange(len(genera))[::-1]
    ylabels, ycolors = [], []
    for k, (g, ph) in enumerate(genera):
        col = S.PHYLUM.get(ph, S.PALETTE["grey"])
        y = ypos[k]
        vals = {f: float(Bm[(Bm["entity"] == g) & (Bm["sub"] == f)]["v"].iloc[0]) for f in fams}
        nrow = Bn[(Bn["entity"] == g) & (Bn["sub"] == "GH1")]
        ntxt = f"n={int(float(nrow['value'].iloc[0]))}" if len(nrow) else ""
        ylabels.append(f"{g}   {ntxt}")
        ycolors.append(col)

        ax.axhline(y, color=S.PALETTE["grid"], lw=0.7, zorder=1)
        ax.plot([min(vals.values()), max(vals.values())], [y, y], color=col,
                lw=1.5, alpha=0.42, solid_capstyle="round", zorder=2)
        for f in fams:
            ax.scatter([vals[f]], [y], marker=MARKER[f], s=27, color=S.FAMILY[f],
                       edgecolors="white", linewidths=0.5, zorder=4)

    ax.set_yticks(ypos)
    ax.set_yticklabels(ylabels, fontsize=S.FS["tick"])
    for lbl, c, (g, ph) in zip(ax.get_yticklabels(), ycolors, genera):
        lbl.set_color(c)
        if ph == "Bacteroidota":
            lbl.set_fontweight("bold")
    ax.set_ylim(-0.9, len(genera) - 0.1)
    ax.set_xscale(yscale, linthresh=ylinthresh)
    ax.set_xticks([0, 1, 2, 5, 10, 20, 50])
    ax.set_xticklabels(["0", "1", "2", "5", "10", "20", "50"])
    ax.set_xlim(-0.9, 62)
    S.style_ax(ax, grid="x", keep_spines=("bottom",))
    ax.set_xlabel("Genus-median GH copy number  (12 genera collapsed from 134 strains; symlog)",
                  fontsize=S.FS["axis"], labelpad=1.5)
    S.panel_label(ax, "D", "Genus level \u2014 the collapse that removes phylum-level significance")
    ax.legend(handles=[Line2D([], [], marker=MARKER[f], linestyle="", markerfacecolor=S.FAMILY[f],
                              markeredgecolor="white", markeredgewidth=0.5, markersize=5.0,
                              label=f"{f} ({S.FAMILY_ROLE[f]})") for f in fams],
              loc="center right", ncol=1, handletextpad=0.4, columnspacing=1.0,
              labelspacing=0.35, bbox_to_anchor=(0.995, 0.45), borderaxespad=0.2)

    # ---------------- 强制科学注记（独立轴，绝不压图）----------------
    st = sd.block("stats")

    def gs_(kind, fam, field):
        r = st[(st["entity"] == kind) & (st["sub"] == fam) & (st["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    D = gs_("phylo_signal", "GH1_presence", "D")
    l2, l3 = gs_("phylo_signal", "GH2", "lambda_ML"), gs_("phylo_signal", "GH3", "lambda_ML")
    pg = {f: gs_("genus_median_phylum_KW", f, "p") for f in fams}
    eps = {f: gs_("strain_phylum_KW", f, "epsilon_squared") for f in fams}

    axN.text(0.0, 1.0,
        "Descriptive pattern + strong phylogenetic signal "
        f"(Fritz\u2013Purvis D = {D:.3f} for GH1 presence/absence; Pagel \u03bb = {l2:.2f} for GH2, {l3:.2f} for GH3).",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", ha="left")
    axN.text(0.0, 0.62,
        "No significance stars are shown by design. Strain-level phylum Kruskal\u2013Wallis is highly significant "
        f"(GH1 \u03b5\u00b2 = {eps['GH1']:.3f}, GH2 \u03b5\u00b2 = {eps['GH2']:.3f}, GH3 \u03b5\u00b2 = {eps['GH3']:.3f}) but that is pseudo-replication within genera;\n"
        f"after collapsing 134 strains to 12 genus medians the same test is not significant (GH1 p = {pg['GH1']:.3f}, GH2 p = {pg['GH2']:.3f}, GH3 p = {pg['GH3']:.3f}).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top", ha="left")
    axN.text(0.0, 0.16,
        "Panels A\u2013C boxes: delivered q1 / median / q3 with min\u2013max whiskers (phylum_family_distribution.tsv); dots: per-strain values (official_strain_cazyme_matrix_round4.tsv). "
        "The symlog axis keeps the true zeros visible rather than dropping them.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top", ha="left")

    S.figure_title(fig,
        "GH1 / GH2 / GH3 copy number across four phyla and twelve genera (134 genomes)",
        x=0.052, y=0.988)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[fig1] PLOT   ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
