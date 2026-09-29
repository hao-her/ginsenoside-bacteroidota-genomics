#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_suppfig1.py — Suppl Fig 1 GH1 阈值曲线（线稿，1200 dpi）
==============================================================
【任务书要求】折线/阶梯图，4 个阈值档，主口径恒为 0；可叠加变体 A/B 作对照线。

【设计思路】
Panel A：四级 i-evalue 阈值（1e-15/1e-10/1e-5/1e-3）下 63 株 Bacteroidales 的 GH1 检出株数。
    主口径（官方全过滤，取自 t5_gh23_control_curves 的 GH1 行）= 0,0,0,0 → 朱红实线贴地；
    变体 A（仅 GH1 家族内 overlap 过滤）= 0,0,3,4 → 灰虚线；
    变体 B（仅 i-evalue、无覆盖度过滤）= 0,0,12,24 → 天蓝点线；
    [1e-3,1) 尾部 25 株以空心灰点作"仅上下文"标记（不参与结论）。
Panel B：阳性对照 GH2（60 株）/ GH3（62 株）在四级阈值下恒定 → 证明管线本身工作正常，
    GH1 的 0 不是技术失败。
Panel C：机制 —— 76 个 loose GH1 域的 覆盖度 × −log10(i-evalue) 散点，画出 cov=0.35 与
    1e-15/1e-5 判据线；说明"放宽阈值后出现的 GH1 信号"是低覆盖度弱亚域（中位 cov 0.169），
    且落在被更强 GH2/GH39 域主导的蛋白上（官方跨家族 overlap 过滤将其归入他家族）。

运行：python plot_suppfig1.py / --plot-only
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import shared_style as S

FIG_NAME = "suppfig1_gh1_threshold_curve"
OUT_DIR = S.SUPP_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 92.0
KIND = "line"                      # 线稿 → 1200 dpi

THRS = ["1e-15", "1e-10", "1e-05", "0.001"]


# ======================================================================
# BUILD
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Suppl Fig 1 — GH1 detection across e-value thresholds (63 Bacteroidales)")
    P = {k: S.rel(S.SRC[k]) for k in ("gh23_curve", "gh1_curve", "frag_inventory")}

    ctl = S.read_upstream("gh23_curve").copy()
    ctl["threshold"] = ctl["threshold"].astype(str)
    for fam in ("GH1", "GH2", "GH3"):
        v = ctl[ctl["family"] == fam].set_index("threshold")["strains_detected"]
        pr = ctl[ctl["family"] == fam].set_index("threshold")["proteins"]
        for i, t in enumerate(THRS):
            if fam == "GH1":   # 主口径只记 GH1 一次
                st.add("A_curve", "official_full_filter", str(i), t, str(int(v.loc[t])), prov=P["gh23_curve"])
            st.add("B_control", fam, str(i), t, str(int(v.loc[t])), prov=P["gh23_curve"])
            st.add("B_control_proteins", fam, str(i), t, str(int(pr.loc[t])), prov=P["gh23_curve"])

    cv = S.read_upstream("gh1_curve").copy()
    cv["threshold"] = cv["threshold"].astype(str)
    for var in ("A_official_full_filters", "B_ievalue_only_no_cov"):
        vv = cv[cv["variant"] == var].set_index("threshold")
        for i, t in enumerate(THRS):
            st.add("A_curve", var, str(i), t, str(int(vv.loc[t, "strains_detected"])), prov=P["gh1_curve"])
        tail = cv[(cv["variant"] == var) & (cv["threshold"].str.contains("tail"))]
        if len(tail):
            st.add("A_tail", var, "", "strains", str(int(tail["strains_detected"].iloc[0])),
                   text=str(tail["threshold"].iloc[0]), prov=P["gh1_curve"])

    fr = S.read_upstream("frag_inventory")
    for idx, (_, r) in enumerate(fr.iterrows()):
        st.add("C_fragments", str(r["accession"]), str(idx), "coverage", str(r["coverage"]),
               text=str(r["target_name"]), prov=P["frag_inventory"])
        st.add("C_fragments", str(r["accession"]), str(idx), "ie", str(r["i_evalue"]),
               text=str(r["target_name"]), prov=P["frag_inventory"])
    st.add("meta", "cov_cutoff", "", "value", "0.35", prov="annotation pipeline §3 (HMM coverage > 0.35)")
    for i, t in enumerate(THRS):
        st.add("meta", "threshold_order", str(i), "name", text=t, prov="task spec")

    st.write(SOURCE_DATA)
    print(f"[supp1] BUILD -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    A = sd.block("A_curve")
    Bc = sd.block("B_control")
    Bp = sd.block("B_control_proteins")
    C = (sd.block("C_fragments")
         .pivot_table(index=["entity", "sub"], columns="field", values="value", aggfunc="first")
         .reset_index())
    C["coverage"] = pd.to_numeric(C["coverage"], errors="coerce")
    C["ie"] = pd.to_numeric(C["ie"], errors="coerce")
    cov_cut = sd.one("meta", "value", entity="cov_cutoff")
    thrs = [t for t in (sd.one_text("meta", "name", entity="threshold_order", sub=str(i))
                        for i in range(4)) if t]
    xs = np.arange(len(thrs))

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    gsA = fig.add_gridspec(1, 1, left=0.055, right=0.360, top=0.860, bottom=0.235)
    axA = fig.add_subplot(gsA[0, 0])
    gsB = fig.add_gridspec(1, 1, left=0.440, right=0.660, top=0.860, bottom=0.235)
    axB = fig.add_subplot(gsB[0, 0])
    gsC = fig.add_gridspec(1, 1, left=0.745, right=0.985, top=0.860, bottom=0.235)
    axC = fig.add_subplot(gsC[0, 0])

    # ---------------- Panel A ----------------
    styles = {"official_full_filter": (S.SEM["not_detected"], "-", 2.0, "o",
                                       "official filters (primary)"),
              "A_official_full_filters": (S.PALETTE["grey"], "--", 1.1, "s",
                                          "GH1-family-only overlap filter"),
              "B_ievalue_only_no_cov": (S.PALETTE["sky"], ":", 1.1, "^",
                                        "i-evalue only, no coverage filter")}
    for var, (col, ls, lw, mk, lab) in styles.items():
        v = A[A["entity"] == var].sort_values("sub", key=lambda s: s.astype(int))
        y = pd.to_numeric(v["value"]).to_numpy()
        axA.plot(xs, y, ls, color=col, lw=lw, marker=mk, ms=4.2, label=lab, zorder=4)
        for xi, yi in zip(xs, y):
            if var == "official_full_filter":
                axA.annotate(f"{yi:g}", (xi, yi), textcoords="offset points", xytext=(0, 7),
                             ha="center", fontsize=S.FS["note"], color=col, fontweight="bold")
    # [1e-3,1) 尾部（仅上下文）
    tail = sd.block("A_tail")
    if len(tail):
        tb = tail[tail["entity"] == "B_ievalue_only_no_cov"]
        if len(tb):
            tv = float(tb["value"].iloc[0])
            axA.scatter([xs[-1] + 0.55], [tv], s=30, facecolors="none",
                        edgecolors=S.PALETTE["grey"], linewidths=0.8, zorder=4)
            axA.annotate(f"[1e-3,1) tail: {tv:g} strains (context only)", (xs[-1] + 0.42, tv),
                         textcoords="offset points", xytext=(-4, -8), fontsize=S.FS["note"],
                         color=S.PALETTE["ink_faint"], va="top", ha="right")
    axA.axhline(0, color=S.SEM["not_detected"], lw=0.7, alpha=0.5, zorder=1)
    axA.set_xticks(xs); axA.set_xticklabels(thrs, fontsize=S.FS["tick"])
    axA.set_xlim(-0.4, len(thrs) + 0.2)
    axA.set_ylim(-1.4, 27.5)
    axA.set_xlabel("GH1 domain i-evalue threshold", fontsize=S.FS["axis"])
    axA.set_ylabel("strains with GH1 detected (of 63)", fontsize=S.FS["axis"])
    axA.legend(fontsize=S.FS["note"], loc="center left", frameon=False)
    S.style_ax(axA, grid="y")
    S.panel_label(axA, "A", "GH1 detection vs threshold")

    # ---------------- Panel B：阳性对照 ----------------
    for fam, col in (("GH2", S.SEM["control"]), ("GH3", S.PALETTE["orange"])):
        v = Bc[Bc["entity"] == fam].sort_values("sub", key=lambda s: s.astype(int))
        pr = Bp[Bp["entity"] == fam].sort_values("sub", key=lambda s: s.astype(int))
        y = pd.to_numeric(v["value"]).to_numpy()
        axB.plot(xs, y, "-", color=col, lw=1.8, marker="o", ms=4.2, zorder=4)
        axB.annotate(f"{fam}: {y[0]:g} strains / {int(float(pr['value'].iloc[0])):g} proteins",
                     (xs[0], y[0]), textcoords="offset points", xytext=(4, -13 if fam == "GH2" else 7),
                     fontsize=S.FS["note"], color=col)
    axB.set_xticks(xs); axB.set_xticklabels(thrs, fontsize=S.FS["tick"])
    axB.set_xlim(-0.4, len(thrs) - 0.6)
    axB.set_ylim(50, 70)
    axB.set_xlabel("i-evalue threshold", fontsize=S.FS["axis"])
    axB.set_ylabel("strains detected (of 63)", fontsize=S.FS["axis"])
    S.style_ax(axB, grid="y")
    S.panel_label(axB, "B", "Positive controls are threshold-stable")

    # ---------------- Panel C：机制（弱亚域）----------------
    neglog = -np.log10(C["ie"].clip(lower=1e-300))
    passes = (C["coverage"] > cov_cut)
    axC.scatter(neglog[~passes], C["coverage"][~passes], s=11, color=S.PALETTE["grey"],
                alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3,
                label="rejected (cov ≤ 0.35)")
    axC.scatter(neglog[passes], C["coverage"][passes], s=16, color=S.SEM["not_detected"],
                edgecolors="white", linewidths=0.4, zorder=4,
                label="cov > 0.35 (survives only at loose thresholds)")
    axC.axhline(cov_cut, color=S.PALETTE["ink_soft"], lw=0.9, ls="--", zorder=2)
    axC.text(0.02, cov_cut + 0.02, f"HMM coverage cutoff = {cov_cut:g}", transform=axC.get_yaxis_transform(),
             fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="bottom")
    for ie, lab in ((1e-15, "1e-15"), (1e-5, "1e-5")):
        axC.axvline(-np.log10(ie), color=S.PALETTE["ink_faint"], lw=0.7, ls=":", zorder=1)
    med_cov = float(C["coverage"].median())
    axC.text(0.98, 0.04, f"n = {len(C)} loose GH1 domains\nmedian coverage = {med_cov:.3f}",
             transform=axC.transAxes, ha="right", va="bottom", fontsize=S.FS["note"],
             color=S.PALETTE["ink_soft"], linespacing=1.35)
    axC.set_xlabel("−log10(i-evalue)", fontsize=S.FS["axis"])
    axC.set_ylabel("HMM coverage", fontsize=S.FS["axis"])
    axC.set_ylim(0, 1.02)
    axC.legend(fontsize=S.FS["note"], loc="upper left", frameon=False)
    S.style_ax(axC, grid="y")
    S.panel_label(axC, "C", "Relaxed-threshold 'GH1' signals are weak sub-domains")

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.055, right=0.985, top=0.150, bottom=0.010)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")
    axN.text(0, 1.0,
        "Under the prescribed dbCAN V5-2 filters GH1 is detected in 0 / 63 Bacteroidales genomes at every threshold from 1e-15 to 1e-3, while GH2 (60 strains) and GH3 (62 strains) "
        "are detected stably — the zero is specific to GH1, not a pipeline failure.",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", fontweight="bold")
    axN.text(0, 0.50,
        "The fragments that appear when filters are relaxed have low HMM coverage (median 0.169) and, at cov > 0.35, sit on proteins dominated by much stronger GH2/GH39 domains "
        "(i-evalue 1e-110 to 1e-214), so the official cross-family overlap filter assigns them to those families.  “Not detected at the prescribed thresholds” ≠ biological absence.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.08,
        "Panel C dotted verticals mark the 1e-15 and 1e-5 i-evalue criteria.  Sources: t5_gh23_control_curves.tsv (official pipeline, panel A primary + panel B), t5_gh1_threshold_curve.tsv (variants A/B + tail), t5_gh1_fragment_inventory.tsv (panel C).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top")

    S.figure_title(fig, "GH1 non-detection is robust to e-value threshold relaxation (63 Bacteroidales genomes)",
                   x=0.055, y=0.988)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[supp1] PLOT  ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
