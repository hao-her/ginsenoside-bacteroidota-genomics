#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_figR2.py — Fig R2 保簇零模型 model C 确认无偏主指标的 PUL 富集
=================================================================
PLOT 阶段：仅从 source_data 三个 TSV 读数。
三面板复合图：
  A  无偏 45 基因组全配对（627 对）森林图：5 指标 × model P/C，
     高斯雾带 = 2000 次置换零分布（ER 尺度），whisker = 95% 置换区间；
     实心 = 通过预注册规则（ER > 2 且 p < 0.01）
  B  GH3 / GH2 / any 指标的置换零分布小提琴（P/C 并排），红菱形 = 观测值
  C  分层（4 层 × core/any）P↔C 哑铃图：两模型 ER 几乎重合
输出：figR2_pul_modelC_unbiased45.{svg,pdf,png(300dpi)}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

import review_style as S

FIG_NAME = "figR2_pul_modelC_unbiased45"
W_MM, H_MM = 183.0, 92.0

METRIC_ORDER = ["core", "GH3", "side", "GH2", "any"]   # 画序：core 最下
METRIC_LABEL = {"any": "Any panel GH", "GH2": "GH2", "side": "Side-chain",
                "GH3": "GH3", "core": "Core (GH1/GH3)"}
LAYER_ORDER = ["L1_Bacteroidales", "L2_Flavobacteriales",
               "L3_Cytophagales", "L3_Sphingobacteriales"]
LAYER_SHORT = {"L1_Bacteroidales": "L1 Bact.",
               "L2_Flavobacteriales": "L2 Flavo.",
               "L3_Cytophagales": "L3 Cyto.",
               "L3_Sphingobacteriales": "L3 Sphingo."}
METRIC_SHORT = {"core": "core", "any": "any"}
MODEL_C = {"P_pairshuffle": S.C["modelP"], "C_ghshift": S.C["modelC"]}
MODEL_L = {"P_pairshuffle": S.C["modelP_l"], "C_ghshift": S.C["modelC_l"]}


def sig(ER: float, p: float) -> bool:
    return (ER > 2.0) and (p < 0.01)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forest", required=True, type=Path)
    ap.add_argument("--perlayer", required=True, type=Path)
    ap.add_argument("--distributions", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    forest = pd.read_csv(a.forest, sep="\t")
    layer = pd.read_csv(a.perlayer, sep="\t")
    dist = pd.read_csv(a.distributions, sep="\t")

    fig = S.new_fig(W_MM, H_MM)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.30, 0.72, 0.98],
                          wspace=0.42, left=0.115, right=0.99, top=0.865, bottom=0.13)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    # 图例（整图顶部，三图共用）
    handles_top = [
        Line2D([], [], marker="o", ls="none", ms=4.6, mfc=S.C["modelP"],
               mec=S.C["modelP"], label="model P (pair permutation)"),
        Line2D([], [], marker="o", ls="none", ms=4.6, mfc=S.C["modelC"],
               mec=S.C["modelC"], label="model C (cluster-preserving shift)"),
        Line2D([], [], marker="o", ls="none", ms=4.6, mfc="white",
               mec=S.C["ink_soft"], label="fails pre-registered rule (ER > 2 & p < 0.01)"),
        Line2D([], [], marker="D", ls="none", ms=4.0, mfc=S.C["vermillion"],
               mec="white", label="observed count (panel B)"),
    ]
    fig.legend(handles=handles_top, loc="upper center", ncol=4, frameon=False,
               handletextpad=0.3, columnspacing=1.0, bbox_to_anchor=(0.53, 0.995))

    # ================= A: forest plot =================
    ypos, ylabels = [], []
    y = 0.0
    for met in METRIC_ORDER:
        y += 1.0
        y_model = {"P_pairshuffle": y + 0.17, "C_ghshift": y - 0.17}
        for model in ("P_pairshuffle", "C_ghshift"):
            r = forest[(forest["metric"] == met) & (forest["model"] == model)].iloc[0]
            ER = r["ER"]
            lo = r["observed"] / r["null_p97_5"]
            hi = r["observed"] / r["null_p02_5"]
            mu = r["observed"] / r["null_mean"]
            sd = r["observed"] * r["null_sd"] / r["null_mean"] ** 2
            xs = np.linspace(lo - 0.55 * sd, hi + 0.55 * sd, 200)
            haze = stats.norm.pdf(xs, mu, sd)
            haze = haze / haze.max() * 0.13
            axA.fill_between(xs, y_model[model] - haze, y_model[model] + haze,
                             color=MODEL_L[model], lw=0, zorder=2)
            axA.plot([lo, hi], [y_model[model]] * 2, color=MODEL_C[model],
                     lw=0.7, zorder=3)
            filled = sig(ER, r["empirical_p"])
            axA.plot(ER, y_model[model], "o", ms=4.6,
                     mfc=MODEL_C[model] if filled else "white",
                     mec=MODEL_C[model], mew=1.0, zorder=4)
            axA.text(hi + 0.10, y_model[model],
                     f"{ER:.2f} (p = {r['empirical_p']:.0e})",
                     fontsize=5.0, color=S.C["ink_soft"], va="center", ha="left")
        ypos.append(y)
        ylabels.append(METRIC_LABEL[met])
        y += 0.52
    axA.axvline(1.0, color=S.C["ink_faint"], lw=0.6, ls="-", zorder=1)
    axA.axvline(2.0, color=S.C["vermillion"], lw=0.8, ls=(0, (4, 2)), zorder=1)
    axA.text(2.02, 0.42, "ER = 2 (pre-registered)", fontsize=5.4,
             color=S.C["vermillion"], ha="left", va="bottom")
    axA.set_yticks(ypos)
    axA.set_yticklabels(ylabels, fontsize=6.3)
    axA.set_xlabel("Enrichment ratio (observed / null mean)")
    axA.set_xlim(0.6, 5.6)
    axA.set_ylim(0.25, y + 0.05)
    axA.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axA)
    S.panel_label(axA, "A", x=-0.34)

    # ================= B: permutation violins =================
    DIST2MODEL = {"P": "P_pairshuffle", "C": "C_ghshift"}
    vmetrics = ["GH3", "GH2", "any"]
    vrows = []
    for met in vmetrics:
        for dmod in ("P", "C"):
            model = DIST2MODEL[dmod]
            sub = dist[(dist["metric"] == met) & (dist["scope"] == "all45")
                       & (dist["model"] == dmod)]
            obs = forest[(forest["metric"] == met) & (forest["model"] == model)].iloc[0]["observed"]
            vrows.append((met, model, sub["null_count"].to_numpy(), obs))
    yy = []
    ycur = 0.0
    pos = []
    for i, (met, model, arr, obs) in enumerate(vrows):
        if i % 2 == 0:
            ycur += 1.0
            base = ycur
        pos.append(base + (0.18 if model == "P_pairshuffle" else -0.18))
        if i % 2 == 1:
            ycur += 0.55
    data = [vr[2] for vr in vrows]
    vp = axB.violinplot(data, positions=pos, orientation="horizontal", widths=0.30,
                        showmeans=False, showmedians=False, showextrema=False)
    for body, (met, model, _, _) in zip(vp["bodies"], vrows):
        body.set_facecolor(MODEL_L[model])
        body.set_edgecolor(MODEL_C[model])
        body.set_linewidth(0.6)
    for pi, (met, model, arr, obs) in zip(pos, vrows):
        q02, q97 = np.percentile(arr, [2.5, 97.5])
        axB.plot([q02, q97], [pi, pi], color=MODEL_C[model], lw=0.8, zorder=3)
        axB.plot(obs, pi, "D", ms=3.4, mfc=S.C["vermillion"], mec="white",
                 mew=0.5, zorder=4)
    met_ticks = [np.mean([pos[2 * i], pos[2 * i + 1]]) for i in range(len(vmetrics))]
    axB.set_yticks(met_ticks)
    axB.set_yticklabels([METRIC_LABEL[m] for m in vmetrics], fontsize=6.3)
    axB.set_xlabel("Null GH-window count (2000 perm.)")
    axB.set_xlim(0, 102)
    axB.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axB)
    S.panel_label(axB, "B", x=-0.30)

    # ================= C: per-layer dumbbell =================
    rows = []
    for lay in LAYER_ORDER:
        for met in ("core", "any"):
            r = {"layer": lay, "metric": met}
            for model in ("P_pairshuffle", "C_ghshift"):
                s = layer[(layer["scope"] == lay) & (layer["metric"] == met)
                          & (layer["model"] == model)].iloc[0]
                r[model + "_ER"] = s["ER"]
                r[model + "_p"] = s["empirical_p"]
            rows.append(r)
    yy = np.arange(len(rows))[::-1].astype(float)
    DODGE = 0.14
    for yi, r in zip(yy, rows):
        pER, cER = r["P_pairshuffle_ER"], r["C_ghshift_ER"]
        axC.plot([pER, cER], [yi, yi], color=S.C["grid"], lw=1.3,
                 solid_capstyle="round", zorder=2)
        for model, ER, yd in (("P_pairshuffle", pER, yi + DODGE),
                              ("C_ghshift", cER, yi - DODGE)):
            filled = sig(ER, r[model + "_p"])
            axC.plot(ER, yd, "o", ms=4.0,
                     mfc=MODEL_C[model] if filled else "white",
                     mec=MODEL_C[model], mew=1.0, zorder=3)
    axC.axvline(2.0, color=S.C["vermillion"], lw=0.8, ls=(0, (4, 2)), zorder=1)
    axC.axvline(1.0, color=S.C["ink_faint"], lw=0.6, zorder=1)
    axC.set_yticks(yy)
    axC.set_yticklabels([f"{LAYER_SHORT[r['layer']]} · {r['metric']}"
                         for r in rows], fontsize=5.6)
    axC.set_xlabel("Enrichment ratio (observed / null mean)")
    axC.set_xlim(0.4, 5.4)
    axC.grid(axis="x", color=S.C["grid"], lw=0.4, zorder=0)
    S.strip_axes(axC)
    S.panel_label(axC, "C", x=-0.56)

    S.save_fig(fig, a.out, FIG_NAME)
    print("saved", FIG_NAME)


if __name__ == "__main__":
    main()
