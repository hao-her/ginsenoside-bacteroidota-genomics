#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_figR1.py — Fig R1 催化对几何预筛（AlphaFold2/ESMFold）
=================================================================
PLOT 阶段：仅从 source_data/figR1_structure_geometry_source_data.tsv 读数。
三面板复合图：
  A  每候选 Cδ–Cδ 距离点图（参照实验值虚线；点形=结构来源）
  B  Cδ–Cδ × min carboxyl O–O 散点（实验参照星形 + 预注册几何判定线）
  C  催化残基 pLDDT 哑铃图（acid/base ↔ nucleophile；判定线 70）
输出：figR1_catalytic_geometry_prescreen.{svg,pdf,png(300dpi)}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import review_style as S

FIG_NAME = "figR1_catalytic_geometry_prescreen"
W_MM, H_MM = 183.0, 80.0


def load(src: Path):
    df = pd.read_csv(src, sep="\t")
    cand = df[df["block"] == "cand"].copy()
    ref = df[df["block"] == "ref"].copy()
    cand["short"] = cand["entity"].map(S.short_wp)
    # GH3 在上（按 rank 升序），GH5 在下
    cand = cand.sort_values(["family", "rank"], ascending=[True, True])
    return cand.reset_index(drop=True), ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-data", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    cand, ref = load(a.source_data)

    fig = S.new_fig(W_MM, H_MM)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 0.92],
                          wspace=0.24, left=0.115, right=0.99, top=0.91, bottom=0.155)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2], sharey=axA)

    fam_c = {"GH3": S.C["GH3"], "GH5": S.C["GH5"]}
    src_mk = {"AlphaFold_DB_v6": "o", "ESMFold_v1_catdomain": "s"}
    y = np.arange(len(cand))[::-1]          # 首行在最上
    n3 = int((cand["family"] == "GH3").sum())

    # ---------- A: per-candidate Cδ–Cδ dot plot ----------
    for (_, r), yi in zip(cand.iterrows(), y):
        axA.plot(r["cdcd_A"], yi, marker=src_mk.get(r["struct_source"], "o"),
                 ms=4.0, mfc=fam_c[r["family"]], mec="white", mew=0.5,
                 ls="none", zorder=3)
    r3 = ref[ref["family"] == "GH3"].iloc[0]
    r5 = ref[ref["family"] == "GH5"].iloc[0]
    axA.axvline(r3["cdcd_A"], color=S.C["GH3"], lw=0.8, ls=(0, (4, 2)), zorder=2)
    axA.axvline(r5["cdcd_A"], color=S.C["GH5"], lw=0.8, ls=(0, (4, 2)), zorder=2)
    axA.text(r3["cdcd_A"], 1.015, "BlBG3 7.39 Å", fontsize=5.6,
             color=S.C["GH3"], ha="center", va="bottom",
             transform=axA.get_xaxis_transform())
    axA.text(r5["cdcd_A"], 1.015, "BcelFp 4.80 Å", fontsize=5.6,
             color=S.C["GH5"], ha="center", va="bottom",
             transform=axA.get_xaxis_transform())
    axA.set_yticks(y)
    axA.set_yticklabels(cand["short"], fontsize=5.5)
    axA.set_xlabel("Catalytic Cδ–Cδ distance (Å)")
    axA.set_xlim(3.6, 9.3)
    axA.set_ylim(-0.8, len(cand) - 0.2)
    axA.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axA)
    # 家族分组（轴右侧色带 + 标签）
    y_gh3_mid = (y[0] + y[n3 - 1]) / 2
    y_gh5_mid = (y[n3] + y[-1]) / 2
    for ymid, ytop, ybot, fam in ((y_gh3_mid, y[0], y[n3 - 1], "GH3"),
                                  (y_gh5_mid, y[n3], y[-1], "GH5")):
        axA.plot([9.0, 9.0], [ybot - 0.4, ytop + 0.4], color=fam_c[fam],
                 lw=2.2, solid_capstyle="butt", clip_on=False, zorder=4)
        axA.text(9.16, ymid, fam, color=fam_c[fam], fontsize=6.5,
                 fontweight="bold", ha="left", va="center")
    S.panel_label(axA, "A", x=-0.46)
    handlesA = [Line2D([], [], marker="o", ls="none", ms=4.0, mfc=S.C["ink_soft"],
                       mec="white", label="AlphaFold DB (full-length)"),
                Line2D([], [], marker="s", ls="none", ms=4.0, mfc=S.C["ink_soft"],
                       mec="white", label="ESMFold (catalytic domain)")]
    axA.legend(handles=handlesA, loc="lower right", frameon=False,
               handletextpad=0.25, borderaxespad=0.15, labelspacing=0.32,
               bbox_to_anchor=(0.995, 0.02))

    # ---------- B: Cδ–Cδ × min O–O scatter ----------
    for fam, sub in cand.groupby("family"):
        pl = np.minimum(sub["plddt_acid_base"], sub["plddt_nucleophile"])
        alpha = 0.45 + 0.55 * (pl - 70) / 30.0
        axB.scatter(sub["cdcd_A"], sub["minOO_A"], s=15, c=[fam_c[fam]],
                    alpha=alpha, edgecolors="white", linewidths=0.5, zorder=3)
    for _, r in ref.iterrows():
        axB.scatter(r["cdcd_A"], r["minOO_A"], marker="*", s=105,
                    c=S.C["ink"], edgecolors="white", linewidths=0.4, zorder=4)
    axB.annotate("BlBG3 (5Z9S)", (r3["cdcd_A"], r3["minOO_A"]),
                 xytext=(5, 6), textcoords="offset points",
                 fontsize=5.8, color=S.C["ink"])
    axB.annotate("BcelFp (6KDD)", (r5["cdcd_A"], r5["minOO_A"]),
                 xytext=(7, -3), textcoords="offset points",
                 fontsize=5.8, color=S.C["ink"])
    axB.axvline(9.0, color=S.C["ink_faint"], lw=0.6, ls=":")
    axB.axhline(8.0, color=S.C["ink_faint"], lw=0.6, ls=":")
    axB.text(8.88, 2.62, "pre-registered cut-offs 9 / 8 Å", fontsize=5.4,
             color=S.C["ink_soft"], ha="right", va="bottom")
    axB.set_xlabel("Catalytic Cδ–Cδ distance (Å)")
    axB.set_ylabel("Min carboxyl O–O distance (Å)")
    axB.set_xlim(3.4, 9.6)
    axB.set_ylim(2.4, 8.4)
    axB.grid(color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axB)
    handles = [Line2D([], [], marker="o", ls="none", ms=4.2, mfc=S.C["GH3"],
                      mec="white", label="GH3 candidate"),
               Line2D([], [], marker="o", ls="none", ms=4.2, mfc=S.C["GH5"],
                      mec="white", label="GH5 candidate"),
               Line2D([], [], marker="*", ls="none", ms=8.5, mfc=S.C["ink"],
                      mec="white", label="Experimental reference")]
    axB.legend(handles=handles, loc="upper left", frameon=False,
               handletextpad=0.25, borderaxespad=0.1, labelspacing=0.32)
    S.panel_label(axB, "B", x=-0.26)

    # ---------- C: pLDDT dumbbell ----------
    for (_, r), yi in zip(cand.iterrows(), y):
        lo, hi = sorted([r["plddt_acid_base"], r["plddt_nucleophile"]])
        axC.plot([lo, hi], [yi, yi], color=S.C["grid"], lw=1.3, zorder=2,
                 solid_capstyle="round")
        axC.plot(r["plddt_acid_base"], yi, "o", ms=3.6, mfc=fam_c[r["family"]],
                 mec="white", mew=0.4, zorder=3)
        axC.plot(r["plddt_nucleophile"], yi, "o", ms=3.6, mfc="white",
                 mec=fam_c[r["family"]], mew=0.9, zorder=3)
    axC.axvline(70, color=S.C["vermillion"], lw=0.8, ls=(0, (4, 2)), zorder=2)
    axC.text(70, 1.015, "cut-off 70", fontsize=5.6, color=S.C["vermillion"],
             ha="left", va="bottom", transform=axC.get_xaxis_transform())
    axC.set_xlabel("Per-residue pLDDT at catalytic Glu")
    axC.set_xlim(67, 101)
    axC.set_ylim(-0.8, len(cand) - 0.2)
    axC.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    plt.setp(axC.get_yticklabels(), visible=False)
    axC.tick_params(axis="y", length=0)
    S.strip_axes(axC)
    handles2 = [Line2D([], [], marker="o", ls="none", ms=3.8, mfc=S.C["ink_soft"],
                       mec="white", label="acid/base Glu"),
                Line2D([], [], marker="o", ls="none", ms=3.8, mfc="white",
                       mec=S.C["ink_soft"], label="nucleophile Glu")]
    axC.legend(handles=handles2, loc="lower left", frameon=False,
               handletextpad=0.25, borderaxespad=0.15, labelspacing=0.32)
    S.panel_label(axC, "C", x=-0.13)

    S.save_fig(fig, a.out, FIG_NAME)
    print("saved", FIG_NAME)


if __name__ == "__main__":
    main()
