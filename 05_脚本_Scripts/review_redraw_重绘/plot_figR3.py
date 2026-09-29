#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_figR3.py — Fig R3 人肠道宏转录组表达证据（473 候选）
=================================================================
PLOT 阶段：仅从 source_data/figR3_metatranscriptome_expression_source_data.tsv 读数。
三面板复合图：
  A  best transcript-fragment identity 的 ECDF（GH3/GH5 分线；90/95% 阈值线；
     注记按表达调用口径：n_ge90_len50 ≥ 1 / n_ge95_len50 ≥ 1 的候选数）
  B  表达调用分类构成（100% 水平堆叠条，段内标 n；GH3 n=417 / GH5 n=56）
  C  best identity × ≥90%-identity fragment count 散点
     （点大小 = ≥95% 片段数；点形 = 家族；标注点名候选）
输出：figR3_metatranscriptome_expression.{svg,pdf,png(300dpi)}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import review_style as S

FIG_NAME = "figR3_metatranscriptome_expression"
W_MM, H_MM = 183.0, 82.0

CALL_ORDER = ["expressed: near-exact transcript (>=95% id, >=80 aa)",
              "expressed: high-identity transcript (>=90% id, >=50 aa)",
              "homolog transcript detected (60-90% id)",
              "no transcript hit (>=60% id)"]
CALL_LABEL = {CALL_ORDER[0]: "near-exact (≥95% id)",
              CALL_ORDER[1]: "high-identity (≥90% id)",
              CALL_ORDER[2]: "homolog (60–90% id)",
              CALL_ORDER[3]: "no hit (≥60% id)"}
CALL_COLOR = {CALL_ORDER[0]: S.C["expr_exact"],
              CALL_ORDER[1]: S.C["expr_high"],
              CALL_ORDER[2]: S.C["expr_homolog"],
              CALL_ORDER[3]: S.C["expr_nohit"]}
FAM_COLOR = {"GH3": S.C["GH3"], "GH5": S.C["GH5"]}

HIGHLIGHT = {
    "WP_195474816.1": ("WP_195474816.1 (B. difficilis, 97.7%)", (-8, -16), "right"),
    "WP_007216976.1": ("WP_007216976.1 (B. cellulosilyticus, 97.5%)", (-9, 4), "right"),
    "WP_005811975.1": ("WP_005811975.1 (B. hominis, 91.8%)", (-8, -18), "right"),
    "WP_010801083.1": ("WP_010801083.1 (tier-1 head,\nhomolog-level only)", (7, 7), "left"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-data", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    df = pd.read_csv(a.source_data, sep="\t")

    fig = S.new_fig(W_MM, H_MM)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.94, 0.86, 1.22],
                        wspace=0.30, left=0.065, right=0.99, top=0.855, bottom=0.16)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    # 顶部统一图例（表达四级分类，B/C 共用色系）
    handles_top = [Patch(facecolor=CALL_COLOR[c], edgecolor="white",
                         label=CALL_LABEL[c]) for c in CALL_ORDER]
    fig.legend(handles=handles_top, loc="upper center", ncol=4, frameon=False,
               handletextpad=0.35, columnspacing=1.0, bbox_to_anchor=(0.53, 0.985))

    # ================= A: ECDF of best identity =================
    for fam in ("GH3", "GH5"):
        v = np.sort(df[df["family"] == fam]["best_pident"].to_numpy())
        y = np.arange(1, len(v) + 1) / len(v)
        axA.step(np.r_[0, v], np.r_[0, y], where="post",
                 color=FAM_COLOR[fam], lw=1.2,
                 label=f"{fam} (n = {len(v)})")
    for x, lab, ha, xx in ((90, "90% id", "right", -1.2), (95, "95% id", "left", 1.2)):
        axA.axvline(x, color=S.C["ink_faint"], lw=0.7, ls=(0, (4, 2)), zorder=1)
    axA.text(88.8, 1.005, "90% id", fontsize=5.6, color=S.C["ink_soft"],
             va="top", ha="right")
    axA.text(96.3, 0.46, "95% id", fontsize=5.6, color=S.C["ink_soft"],
             va="center", ha="left")
    # 表达调用口径计数（与 Table S9 / 正文一致）：
    # ≥90% = near-exact + high-identity；≥95% = near-exact
    n90 = int(df["metatx_expression_call"].isin(CALL_ORDER[:2]).sum())
    n95 = int((df["metatx_expression_call"] == CALL_ORDER[0]).sum())
    ntot = len(df)
    axA.text(3, 0.86,
             f"with ≥90%-id fragments: {n90}/{ntot} ({n90 / ntot:.0%})\n"
             f"with ≥95%-id fragments: {n95}/{ntot} ({n95 / ntot:.0%})",
             fontsize=5.8, color=S.C["ink"], va="top", ha="left")
    axA.set_xlabel("Best transcript-fragment identity (%)")
    axA.set_ylabel("Cumulative fraction of candidates")
    axA.set_xlim(-2, 105)
    axA.set_ylim(0, 1.03)
    axA.grid(color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axA)
    axA.legend(loc="lower right", frameon=False, handletextpad=0.4,
               borderaxespad=0.1)
    S.panel_label(axA, "A", x=-0.24)

    # ================= B: composition 100% stacked bars =================
    ylabels, ypos = [], []
    y = 0.0
    for fam in ("GH3", "GH5"):
        sub = df[df["family"] == fam]
        counts = np.array([int((sub["metatx_expression_call"] == c).sum())
                           for c in CALL_ORDER], dtype=float)
        fracs = counts / counts.sum()
        y += 1.0
        left = 0.0
        for frac, cnt, call in zip(fracs, counts, CALL_ORDER):
            if frac > 0:
                axB.barh(y, frac, left=left, height=0.52, color=CALL_COLOR[call],
                         edgecolor="white", lw=0.5, zorder=3)
                if frac >= 0.055:
                    axB.text(left + frac / 2, y, str(int(cnt)), ha="center",
                             va="center", fontsize=5.6, fontweight="bold", zorder=4,
                             color="white" if call in (CALL_ORDER[0], CALL_ORDER[2])
                             else S.C["ink"])
            left += frac
        axB.text(1.02, y, f"n = {len(sub)}", fontsize=5.6,
                 color=S.C["ink_soft"], va="center", ha="left", zorder=4,
                 transform=axB.get_yaxis_transform())
        ypos.append(y)
        ylabels.append(fam)
        y += 0.85
    axB.set_yticks(ypos)
    axB.set_yticklabels(ylabels, fontsize=6.5)
    axB.set_xlabel("Fraction of candidate passers")
    axB.set_xlim(0, 1.0)
    axB.set_ylim(0.2, y)
    axB.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axB)
    S.panel_label(axB, "B", x=-0.14)

    # ================= C: scatter identity × fragment count =================
    rng = np.random.default_rng(20260927)
    for call in CALL_ORDER:
        sub = df[df["metatx_expression_call"] == call]
        if len(sub) == 0:
            continue
        x = sub["best_pident"].to_numpy()
        yv = sub["n_ge90_len50"].to_numpy(float)
        yj = yv + rng.uniform(-0.22, 0.22, len(sub))
        sizes = 3.5 + sub["n_ge95_len50"].to_numpy(float) * 1.1
        fams = sub["family"].to_numpy()
        for fam in ("GH3", "GH5"):
            m = fams == fam
            if m.any():
                axC.scatter(x[m], yj[m], s=sizes[m], c=[CALL_COLOR[call]],
                            alpha=0.55 if fam == "GH3" else 0.8,
                            edgecolors="none", zorder=3,
                            marker="o" if fam == "GH3" else "^")
    for wp, (lab, (dx, dy), ha) in HIGHLIGHT.items():
        r = df[df["protein_id"] == wp]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        axC.scatter([r["best_pident"]], [r["n_ge90_len50"]], s=36,
                    facecolors="none", edgecolors=S.C["ink"], linewidths=0.9,
                    zorder=5)
        axC.annotate(lab, xy=(r["best_pident"], r["n_ge90_len50"]),
                     xytext=(dx, dy), textcoords="offset points",
                     fontsize=5.2, ha=ha, color=S.C["ink"],
                     arrowprops=dict(arrowstyle="-", lw=0.5,
                                     color=S.C["ink_faint"]))
    for x in (90, 95):
        axC.axvline(x, color=S.C["ink_faint"], lw=0.7, ls=(0, (4, 2)), zorder=1)
    axC.set_xlabel("Best transcript-fragment identity (%)")
    axC.set_ylabel("# ≥90%-id frags.\n(size: # ≥95%-id)")
    axC.set_xlim(-4, 108)
    axC.set_ylim(-1.8, 36)
    axC.grid(color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axC)
    handlesC = [Line2D([], [], marker="o", ls="none", ms=4.0, mfc=S.C["ink_soft"],
                       mec="none", label="GH3"),
                Line2D([], [], marker="^", ls="none", ms=4.4, mfc=S.C["ink_soft"],
                       mec="none", label="GH5")]
    axC.legend(handles=handlesC, loc="upper left", frameon=False,
               handletextpad=0.3, borderaxespad=0.1)
    S.panel_label(axC, "C", x=-0.16)

    S.save_fig(fig, a.out, FIG_NAME)
    print("saved", FIG_NAME)


if __name__ == "__main__":
    main()
