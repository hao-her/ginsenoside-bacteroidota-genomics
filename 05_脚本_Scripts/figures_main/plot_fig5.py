#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig5.py — Fig 5 候选酶筛选漏斗 + 同源性散点
==================================================
【任务书要求】漏斗图（484 → 417 → 473 排序 → top20 → top5）+
同源性散点图（横轴=与参照最高同一性%，颜色=属）。

【设计思路】
Panel A：阶梯漏斗。宽度 ∝ sqrt(计数)（保证 5 这种小阶段仍可见，且不扭曲相对大小的可读性）。
    阶段：GH3 池 484 → GH3 催化残基完整 417（352 complete + 65 pairwise-rescued）
          → 并入 GH5 56 后按预注册 C6 排序的过筛集 473 → tier-1（rank ≤ 20）→ 报告点名的 5 个湿实验候选。
    侧枝：GH5 池 68 → 56；GH66 池 11 → 催化残基 NOT COMPUTABLE（PStDex 序列未公开，斜纹留空，红线 R1）。
Panel B：同源性散点。x = 与最佳参照酶的同一性(%)；y = 属（分类行，按候选数降序）——
    颜色与行位置冗余编码属，故即使两种色在某色觉下接近也不丢信息；
    实心 = 位于 SusC–SusD（PUL）窗内，空心 = 窗外；圆 = GH3，菱 = GH5。
    点名 4 个关键候选（含 GH5_25 那条）并画 S4a 上限参考线（无 ≥98% 命中）。
Panel C：催化判定构成（complete / pairwise-rescued / incomplete / mapping-failure / not-computable）。

红线：不对候选做"能催化 Rb1→CK"的功能断言；仅呈报事实性筛选属性。

运行：python plot_fig5.py / --plot-only
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import textwrap
import numpy as np
import pandas as pd

import shared_style as S

FIG_NAME = "fig5_candidate_screening_funnel"
OUT_DIR = S.MAIN_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 150.0
KIND = "main"

# 报告 §1.4 点名的关键候选（事实性描述，非功能断言）
HIGHLIGHT = ["WP_195474816.1", "WP_005811975.1", "WP_258914935.1", "WP_010801083.1"]


# ======================================================================
# BUILD
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Fig 5 — single-enzyme candidate screening funnel and homology")
    P = {k: S.rel(S.SRC[k]) for k in ("gh3_cand", "gh5_cand", "gh66_inv", "ranked", "shortlist")}

    g3 = S.read_upstream("gh3_cand")
    g5 = S.read_upstream("gh5_cand")
    g66 = S.read_upstream("gh66_inv")
    rk = S.read_upstream("ranked")
    sl = S.read_upstream("shortlist")

    # --- A_funnel：各阶段计数（全部为交付表的行计数/字段计数，不重算统计）---
    def status_counts(df):
        return df["cat_status"].value_counts().to_dict()

    s3, s5 = status_counts(g3), status_counts(g5)
    gh3_pass = s3.get("complete", 0) + s3.get("complete_pairwise_rescued", 0)
    gh5_pass = s5.get("complete", 0) + s5.get("complete_pairwise_rescued", 0)

    stages = [
        ("gh3_pool", len(g3), "GH3 pool (63 Bacteroidales genomes)", P["gh3_cand"]),
        ("gh3_s1_pass", gh3_pass, "GH3 catalytic residues complete (S1)", P["gh3_cand"]),
        ("ranked_passers", len(rk), "merged GH3+GH5 passers, C6-ranked", P["ranked"]),
        ("tier1_top20", int((rk["tier"] == "tier1_top20").sum()), "tier-1 (C6 rank ≤ 20)", P["ranked"]),
        ("highlighted", 5, "report-highlighted wet-lab candidates", P["ranked"]),
    ]
    for name, n, desc, prov in stages:
        st.add("A_funnel", name, "", "n", str(n), text=desc, prov=prov)
    st.add("A_funnel", "gh3_pool", "", "n_strains", "63", prov=P["gh3_cand"])
    st.add("A_funnel", "gh5_pool", "", "n", str(len(g5)), text="GH5 pool — 31 Bacteroidota genomes", prov=P["gh5_cand"])
    st.add("A_funnel", "gh5_s1_pass", "", "n", str(gh5_pass), text="GH5 catalytic residues complete (S1)", prov=P["gh5_cand"])
    st.add("A_funnel", "gh66_pool", "", "n", str(len(g66)), text="GH66 pool — 31 Bacteroidota genomes", prov=P["gh66_inv"])
    st.add("A_funnel", "gh66_not_computable", "", "n", str(len(g66)),
           text="GH66 catalytic residues NOT COMPUTABLE (PStDex sequence not public; no curated ACT_SITE)", prov=P["gh66_inv"])
    for k, v in s3.items():
        st.add("C_verdict", "GH3", "", k, str(v), prov=P["gh3_cand"])
    for k, v in s5.items():
        st.add("C_verdict", "GH5", "", k, str(v), prov=P["gh5_cand"])
    st.add("C_verdict", "GH66", "", "not_computable", str(len(g66)), prov=P["gh66_inv"])
    st.add("A_funnel", "shortlist_clusters", "", "n", str(len(sl)),
           text="non-redundant (≥95% identity) wet-lab shortlist clusters", prov=P["shortlist"])

    # --- B_scatter：473 过筛候选的同源性 / 属 / PUL / 家族 ----------------
    for _, r in rk.iterrows():
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "identity", str(r["identity_best_ref"]),
               prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "tlen", str(r["tlen"]), prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "family", str(r["family"]), prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "pul", str(r["pul_window"]), prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "s4a_max", str(r["s4a_max_identity"]), prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "s4a_reported", str(r["s4a_reported"]), prov=P["ranked"])
        st.add("B_scatter", str(r["protein_id"]), str(r["genera"]), "rank", str(r["rank"]), prov=P["ranked"])
    st.add("meta", "s4a_threshold", "", "value", "98.0",
           text="S4a: a candidate ≥98% identical to a reported enzyme would be excluded as 'already reported'", prov="round9 criteria C6/S4a")
    for i, h in enumerate(HIGHLIGHT):
        st.add("meta", "highlight", str(i), "name", text=h, prov="report_round9 §1.4 / handoff §5")

    st.write(SOURCE_DATA)
    print(f"[fig5] BUILD  -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    F = sd.block("A_funnel")
    B = (sd.block("B_scatter")
         .pivot_table(index=["entity", "sub"], columns="field", values="value", aggfunc="first")
         .reset_index())
    C = sd.block("C_verdict")
    B["identity"] = pd.to_numeric(B["identity"], errors="coerce")

    def fn(name, field="n"):
        r = F[(F["entity"] == name) & (F["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    gsA = fig.add_gridspec(1, 1, left=0.028, right=0.415, top=0.900, bottom=0.160)
    axA = fig.add_subplot(gsA[0, 0])
    gsB = fig.add_gridspec(1, 1, left=0.560, right=0.985, top=0.900, bottom=0.430)
    axB = fig.add_subplot(gsB[0, 0])
    gsC = fig.add_gridspec(1, 1, left=0.560, right=0.985, top=0.300, bottom=0.185)
    axC = fig.add_subplot(gsC[0, 0])

    # ---------------- Panel A：阶梯漏斗（描述文字收在本轴内）----------------
    order = ["gh3_pool", "gh3_s1_pass", "ranked_passers", "tier1_top20", "highlighted"]
    counts = [fn(o) for o in order]
    wmax = np.sqrt(max(counts))
    yy = np.arange(len(order))[::-1] * 1.0
    half = [np.sqrt(c) / wmax * 0.42 for c in counts]
    cols = [S.SEM["pass"]] * 4 + [S.SEM["highlight"]]

    for i, (o, c, hw, col) in enumerate(zip(order, counts, half, cols)):
        y = yy[i]
        axA.add_patch(plt.Rectangle((0.5 - hw, y - 0.30), 2 * hw, 0.60, facecolor=col,
                                    edgecolor="white", linewidth=0.6, alpha=0.92, zorder=3))
        if i < len(order) - 1:
            hw2 = half[i + 1]
            axA.add_patch(plt.Polygon([[0.5 - hw, y - 0.30], [0.5 + hw, y - 0.30],
                                       [0.5 + hw2, yy[i + 1] + 0.30], [0.5 - hw2, yy[i + 1] + 0.30]],
                                      facecolor=col, alpha=0.26, edgecolor="none", zorder=2))
        desc = F[(F["entity"] == o) & (F["field"] == "n")]["text"].iloc[0]
        axA.text(0.5, y, f"{c:g}", ha="center", va="center", fontsize=S.FS["axis"],
                 color="white" if hw > 0.09 else S.PALETTE["ink"], fontweight="bold", zorder=5)
        axA.text(0.96, y, textwrap.fill(desc, 34), ha="left", va="center",
                 fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], zorder=5, linespacing=1.25)

    # 侧枝（GH5 并入 / GH66 不可算 / 非冗余簇）—— 全部在本轴右侧栏内
    axA.text(0.04, -0.95,
             textwrap.fill(f"side track:  GH5 pool {fn('gh5_pool'):g} \u2192 {fn('gh5_s1_pass'):g} pass, merged into the ranked set.  "
                           f"GH66 pool {fn('gh66_pool'):g} \u2192 catalytic residues NOT COMPUTABLE (PStDex sequence not public) \u2014 excluded, not imputed.", 74),
             fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top", ha="left",
             linespacing=1.35)
    axA.text(0.04, -1.95,
             textwrap.fill(f"non-redundant (\u226595% identity) shortlist: {fn('shortlist_clusters'):g} clusters", 74),
             fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top", ha="left")

    axA.set_xlim(0, 2.35)
    axA.set_ylim(-2.9, len(order) - 0.2)
    axA.set_xticks([]); axA.set_yticks([])
    for sp in ("left", "bottom", "right", "top"):
        axA.spines[sp].set_visible(False)
    S.panel_label(axA, "A", "Screening funnel (sqrt-count width scale)")

    # ---------------- Panel B：同源性散点 ----------------
    genera = B["sub"].value_counts()
    gl = list(genera.index)
    ypos = {g: (len(gl) - 1 - i) for i, g in enumerate(gl)}
    rng = np.random.default_rng(7)
    for g in gl:
        sub = B[B["sub"] == g]
        y = ypos[g] + rng.uniform(-0.30, 0.30, len(sub))
        col = S.GENUS.get(g, S.PALETTE["grey"])
        for (xi, yi, fam, pul) in zip(sub["identity"], y, sub["family"], sub["pul"]):
            axB.scatter([xi], [yi], marker="o" if fam == "GH3" else "D", s=13,
                        facecolors=col if pul == "in" else "none",
                        edgecolors=col, linewidths=0.5, zorder=3, alpha=0.85)

    # 点名候选：只画环，说明集中放在右下角文字块
    hl_lines = []
    for i in range(4):
        h = sd.one_text("meta", "name", entity="highlight", sub=str(i))
        if not h:
            continue
        r = B[B["entity"] == h]
        if not len(r):
            continue
        xi = float(r["identity"].iloc[0]); y0 = ypos[r["sub"].iloc[0]]
        axB.scatter([xi], [y0], s=110, facecolors="none", edgecolors=S.SEM["highlight"],
                    linewidths=1.2, zorder=6)
        hl_lines.append(f"{h}  {xi:.2f}%  ({r['sub'].iloc[0]}"
                        + (", GH5_25" if r["family"].iloc[0] == "GH5" else "") + ")")

    s4a_max = float(B["s4a_max"].max())
    thr = sd.one("meta", "value", entity="s4a_threshold")
    n_rep = int((B["s4a_reported"].str.lower() == "true").sum())
    axB.axvline(s4a_max, color=S.PALETTE["ink_faint"], lw=0.8, ls=(0, (3, 2)), zorder=2)

    axB.set_yticks(list(ypos.values()))
    axB.set_yticklabels([f"{g}  (n={genera[g]})" for g in gl], fontsize=S.FS["tick"])
    for lbl, g in zip(axB.get_yticklabels(), gl):
        lbl.set_color(S.GENUS.get(g, S.PALETTE["ink_soft"]))
        if g in ("Bacteroides", "Parabacteroides"):
            lbl.set_fontweight("bold")
    axB.set_ylim(-0.8, len(gl) - 0.2)
    axB.set_xlim(8, 72)
    axB.set_xticks([10, 20, 30, 40, 50, 60, 70])
    axB.set_xlabel("highest identity to a reference enzyme (%)", fontsize=S.FS["axis"])
    S.style_ax(axB, grid="x", keep_spines=("bottom", "left"))
    S.panel_label(axB, "B", "Homology of the 473 ranked passers  (colour = genus · fill = PUL window · shape = family)")

    # ---------------- Panel C：催化判定构成 ----------------
    cats = ["complete", "complete_pairwise_rescued", "incomplete", "mapping_failure", "not_computable"]
    cat_col = {"complete": S.SEM["pass"], "complete_pairwise_rescued": S.SEM["rescued"],
               "incomplete": S.PALETTE["grey"], "mapping_failure": S.PALETTE["grey_l"],
               "not_computable": S.PALETTE["grey_l"]}
    fams = ["GH3", "GH5", "GH66"]
    yc = np.arange(len(fams))[::-1]
    left = np.zeros(len(fams))
    for c in cats:
        vals = np.array([float(C[(C["entity"] == f) & (C["field"] == c)]["value"].iloc[0])
                         if len(C[(C["entity"] == f) & (C["field"] == c)]) else 0.0 for f in fams])
        axC.barh(yc, vals, left=left, height=0.55, color=cat_col[c],
                 edgecolor="white", linewidth=0.5, zorder=3,
                 hatch="////" if c == "not_computable" else None)
        for i, v in enumerate(vals):
            if v > 25 and c not in ("mapping_failure", "not_computable"):
                axC.text(left[i] + v / 2, yc[i], f"{v:g}", ha="center", va="center",
                         fontsize=S.FS["note"], color="white", fontweight="bold")
        left += vals
    axC.set_xlim(0, 980)
    axC.set_yticks(yc)
    axC.set_yticklabels(fams, fontsize=S.FS["tick"])
    axC.set_xlabel("candidates by catalytic-residue verdict", fontsize=S.FS["axis"])
    S.style_ax(axC, grid="x", keep_spines=("bottom",))
    S.panel_label(axC, "C", "Catalytic-residue verdict")
    axC.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=cat_col[c], ec="white",
                                      hatch="////" if c == "not_computable" else None,
                                      label=c.replace("_", " ")) for c in cats],
               loc="center left", bbox_to_anchor=(0.55, 0.5), fontsize=S.FS["note"],
               ncol=1, columnspacing=0.9, handletextpad=0.3, labelspacing=0.30)

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.028, right=0.985, top=0.130, bottom=0.006)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")
    axN.text(0, 1.0,
        "Candidates are reported only as: catalytic residues complete (homology-mapped, structure-anchored to PDB 5Z9S / 6KDD) + family/subfamily assignment + length + "
        "culturability + literature status + PUL co-localisation.  NO claim is made that any candidate catalyses Rb1\u2192CK or any step of it; functional verification requires wet-lab work.",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", fontweight="bold")
    axN.text(0, 0.46,
        f"Panel B: max identity to any reported enzyme = {s4a_max:.2f}% ({n_rep} candidates \u2265 {thr:g}% \u2192 none is \u201calready reported\u201d).  "
        + "  ".join(hl_lines) + ".",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.20,
        "Parabacteroides is a literature-blank genus (no saponin-conversion report), giving it the highest wet-lab novelty; the two highest-identity Bacteroides candidates (66.75% / 66.05% to BglNg-767) "
        "rank 26/28 under the pre-registered C6 order because that order prioritises literature-blank taxa \u2014 both orderings are shown (funnel = C6 rank; panel B = identity).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")

    S.figure_title(fig, "From 484 GH3 proteins to 5 wet-lab candidates: catalytic-residue screening and homology landscape",
                   x=0.030, y=0.988)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[fig5] PLOT   ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
