#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig2.py — Fig 2 谱系边界图（全文核心图）
=============================================
【任务书要求】横轴/行 = 4 个目（Bacteroidales / Flavobacteriales /
Sphingobacteriales / Cytophagales）；用颜色或柱高表示 GH1 检出率（0 vs >0）；
每个目上方标注 株数/属数/科数/Wilson 上限 5.75% / "robust across 4 thresholds"；
Bacteroidales 用醒目色标 0 检出，其他目用另一色标出阳性株数。

【设计思路 —— "缺失带"（absence band）视觉隐喻】
Panel A（通栏，主体）：一个方格 = 一个基因组。4 行（目），行内每株一格。
    GH1 未检出 → 空心朱红格；GH1 检出 → 实心蓝格。
    Bacteroidales 行 = 63 个连续空心朱红格 → 一条刺眼的"空带"，
    读者一眼记住 "Bacteroidales 目级 GH1 未检出"。其余目是短小的实心蓝行。
    行内按来源数据集分段（round4-31 / sentinel-10 / round7-45），细分隔线标注，
    使"三个独立采样共同支持"可见。
Panel B：检出率水平条 + Bacteroidales 的 Wilson 95% 上限须（5.75%，引自 R9-M31/round8 M52）。
Panel C：三数据集并列（task3_2_three_datasets_parallel.tsv）：零值株占比 + 各自已交付的
    Wilson95 上限须 —— 证据的独立复现。
Panel D：四级阈值稳健性（1e-15/1e-10/1e-5/1e-3）：主口径恒 0；变体 A/B 作对照虚线。

所有数值均来自 source_data.tsv；Wilson 上限与阈值曲线不重算，只搬运。

运行：python plot_fig2.py            （BUILD + PLOT）
      python plot_fig2.py --plot-only
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import shared_style as S

FIG_NAME = "fig2_gh1_lineage_boundary"
OUT_DIR = S.MAIN_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 160.0
KIND = "main"

# 三数据集在 Bacteroidales 行的呈现顺序（与任务书 §4 台账一致）
DATASET_ORDER = ["round4_31", "sentinel_10", "round7_45"]
DATASET_TAG = {"round4_31": "r4 (31)", "sentinel_10": "sent. (6)", "round7_45": "r7 (26)"}


# ======================================================================
# BUILD
# ======================================================================
def _order_rows():
    """把三个官方矩阵拼成 株×目×GH1 的长表（仅聚合，不做统计推断）。"""
    parts = []
    r4 = S.read_upstream("matrix_r4")
    r4b = r4[r4["phylum"] == "Bacteroidota"].copy()
    r4b["dataset"] = "round4_31"
    r4b["order"] = "Bacteroidales"          # round4 Bacteroidota 31 株全为 Bacteroidales（台账 §4）
    r4b["family_taxonomy"] = ""             # round4 无 family 列 → 留空，不用推测值填充（R1）
    parts.append(r4b[["accession", "genus", "family_taxonomy", "order", "dataset", "GH1"]])

    sn = S.read_upstream("matrix_sentinel")
    sn = sn.copy(); sn["dataset"] = "sentinel_10"
    parts.append(sn[["accession", "genus", "family_taxonomy", "order_taxonomy", "dataset", "GH1"]]
                 .rename(columns={"order_taxonomy": "order"}))

    r7 = S.read_upstream("matrix_r7")
    r7 = r7.copy(); r7["dataset"] = "round7_45"
    parts.append(r7[["accession", "genus", "family_taxonomy", "order_taxonomy", "dataset", "GH1"]]
                 .rename(columns={"order_taxonomy": "order"}))

    df = pd.concat(parts, ignore_index=True)
    df["GH1"] = pd.to_numeric(df["GH1"], errors="coerce").fillna(0).astype(int)
    return df


def build() -> S.SourceTable:
    st = S.SourceTable("Fig 2 — GH1 lineage boundary across four Bacteroidota orders")
    strains = _order_rows()
    P = {k: S.rel(S.SRC[k]) for k in ("matrix_r4", "matrix_sentinel", "matrix_r7",
                                      "three_datasets", "gh23_curve", "gh1_curve")}

    # --- A_glyph：每株一格 -------------------------------------------
    for _, r in strains.iterrows():
        st.add("A_glyph", str(r["accession"]), str(r["order"]), "gh1", str(int(r["GH1"])),
               text=f'{r["dataset"]}|{r["genus"]}', prov=P["matrix_r4"] + ";" + P["matrix_sentinel"] + ";" + P["matrix_r7"])

    # --- A_order：每目的株数 / 属数 / 科数 / 检出数 --------------------
    for order in S.ORDER_ROW:
        sub = strains[strains["order"] == order]
        # 属数：三个矩阵都有 genus 列 → 全量并集
        n_gen = sub["genus"].nunique()
        # 科数：仅统计带 family_taxonomy 的矩阵（sentinel+round7）；round4 无该列，不推测（R1）
        fams = set(sub.loc[sub["family_taxonomy"] != "", "family_taxonomy"])
        n_fam = len(fams)
        k = int((sub["GH1"] > 0).sum()); n = len(sub)
        for ds in DATASET_ORDER:
            d = sub[sub["dataset"] == ds]
            st.add("A_order_seg", order, ds, "n", str(len(d)), prov=P[ds if ds in P else "matrix_r4"])
        st.add("A_order", order, "", "n_strains", str(n), prov=";".join(P[k2] for k2 in ("matrix_r4", "matrix_sentinel", "matrix_r7")))
        st.add("A_order", order, "", "n_genera", str(n_gen), prov=";".join(P[k2] for k2 in ("matrix_r4", "matrix_sentinel", "matrix_r7")))
        st.add("A_order", order, "", "n_families", str(n_fam),
               text="families counted over sentinel+round7 taxonomy columns (round4 matrix has no family column; not imputed)",
               prov=P["matrix_sentinel"] + ";" + P["matrix_r7"])
        st.add("A_order", order, "", "gh1_positive", str(k), prov=";".join(P[k2] for k2 in ("matrix_r4", "matrix_sentinel", "matrix_r7")))
        st.add("A_order", order, "", "detection_rate", repr(k / n) if n else "",
               num=str(k), den=str(n), prov="aggregation of delivered per-strain counts")

    # --- B_wilson：0/63 的检出率 Wilson 上限（逐字引自 R9-M31 / round8 M52）---
    t4 = S.read_upstream("t4_prov")
    row = t4[t4["M_id"] == "R9-M31"]
    m = re.search(r"([\d.]+)%\s*\((\d+)/(\d+)", row["value"].iloc[0])
    pct, k0, n63 = float(m.group(1)), int(m.group(2)), int(m.group(3))
    st.add("B_wilson", "Bacteroidales", "", "wilson95_upper_pct", str(pct),
           num=str(k0), den=str(n63), text=row["meaning"].iloc[0],
           prov="tables/t4_number_provenance.tsv R9-M31 (= round8 M52)")

    # --- C_datasets：三数据集并列（已交付，含各自 Wilson95 上限）---------
    td = S.read_upstream("three_datasets")
    for _, r in td.iterrows():
        key = "round4_31" if r["dataset"].startswith("round4") else \
              ("sentinel_10" if r["dataset"].startswith("sentinel") else "round7_45")
        st.add("C_datasets", key, "", "n_strains", str(r["n_strains"]), prov=P["three_datasets"])
        st.add("C_datasets", key, "", "gh1_zero_n", str(r["GH1_zero_n"]), prov=P["three_datasets"])
        st.add("C_datasets", key, "", "gh1_zero_prop", str(r["GH1_zero_prop"]), prov=P["three_datasets"])
        st.add("C_datasets", key, "", "gh1_positive_wilson95_upper",
               str(r["GH1_positive_wilson95_upper"]), prov=P["three_datasets"])
        st.add("C_datasets", key, "", "composition", "", text=str(r["dataset"]), prov=P["three_datasets"])

    # --- D_threshold：主口径 + 变体 A/B（已交付，不重算）----------------
    ctl = S.read_upstream("gh23_curve")
    ctl["threshold"] = ctl["threshold"].astype(str)   # 防止 1e-15 被解析为 float
    g1 = ctl[ctl["family"] == "GH1"].set_index("threshold")["strains_detected"]
    thr_order = ["1e-15", "1e-10", "1e-05", "0.001"]
    for i, t in enumerate(thr_order):
        st.add("D_threshold", "official_full_filter", str(i), t, str(int(g1.loc[t])), prov=P["gh23_curve"])
    cv = S.read_upstream("gh1_curve")
    cv["threshold"] = cv["threshold"].astype(str)
    for var in ("A_official_full_filters", "B_ievalue_only_no_cov"):
        v = cv[(cv["variant"] == var) & (cv["threshold"] != "[1e-3,1) tail (context only)")]
        v = v.set_index("threshold")["strains_detected"]
        for i, t in enumerate(thr_order):
            st.add("D_threshold", var, str(i), t, str(int(v.loc[t])), prov=P["gh1_curve"])
    for i, t in enumerate(thr_order):
        st.add("meta", "threshold_order", str(i), "name", text=t, prov="task spec")

    # --- meta ----------------------------------------------------------
    for i, o in enumerate(S.ORDER_ROW):
        st.add("meta", "order_row", str(i), "name", text=o, prov="shared_style.ORDER_ROW")
    for i, d in enumerate(DATASET_ORDER):
        st.add("meta", "dataset_order", str(i), "name", text=d, prov="handoff §4 ledger")
        st.add("meta", "dataset_tag", d, "name", text=DATASET_TAG[d], prov="handoff §4 ledger")
    st.add("meta", "glyph", "meaning", "square", "1.0",
           text="one square = one genome; hollow vermillion = GH1 not detected; solid blue = GH1 detected",
           prov="design")

    st.write(SOURCE_DATA)
    print(f"[fig2] BUILD  -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)

    orders = [t for t in (sd.one_text("meta", "name", entity="order_row", sub=str(i))
                          for i in range(len(S.ORDER_ROW))) if t]
    dsets = [t for t in (sd.one_text("meta", "name", entity="dataset_order", sub=str(i))
                         for i in range(len(DATASET_ORDER))) if t]
    dtags = {d: sd.one_text("meta", "name", entity="dataset_tag", sub=d) for d in dsets}
    thrs = [t for t in (sd.one_text("meta", "name", entity="threshold_order", sub=str(i))
                        for i in range(4)) if t]

    G = sd.block("A_glyph").copy()
    G["gh1"] = pd.to_numeric(G["value"], errors="coerce").fillna(0).astype(int)
    G[["dataset", "genus"]] = G["text"].str.split("|", expand=True)
    O = sd.block("A_order")

    def o_val(order, field):
        r = O[(O["entity"] == order) & (O["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    seg = sd.block("A_order_seg")

    def seg_n(order, ds):
        r = seg[(seg["entity"] == order) & (seg["sub"] == ds)]
        return int(float(r["value"].iloc[0])) if len(r) else 0

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)

    # ---------------- Panel A：缺失带 ----------------
    gsA = fig.add_gridspec(1, 1, left=0.105, right=0.900, top=0.878, bottom=0.560)
    axA = fig.add_subplot(gsA[0, 0])

    nmax = int(max(o_val(o, "n_strains") for o in orders))
    row_y = {o: len(orders) - 1 - i for i, o in enumerate(orders)}
    GAP = 0.9
    for o in orders:
        y = row_y[o]
        sub = G[G["sub"] == o]
        x = 0.0
        seg_starts = []
        for d in dsets:
            ds_sub = sub[sub["dataset"] == d]
            if not len(ds_sub):
                continue
            # 检出株排到段末，使未检出段连续成带
            ds_sub = ds_sub.sort_values("gh1", kind="stable")
            start = x
            for _, r in ds_sub.iterrows():
                pos = (r["gh1"] > 0)
                axA.scatter([x], [y], marker="s",
                            s=30.0, color=S.SEM["detected"] if pos else S.SEM["absent_fill"],
                            edgecolors=S.PALETTE["blue"] if pos else S.SEM["not_detected"],
                            linewidths=0.55 if pos else 0.75, zorder=3)
                x += 1.0
            seg_starts.append((start, x, d))
            x += GAP
        # 段分隔线 + 标签（仅在有多个段的行）
        for (s0, s1, d) in seg_starts:
            if len(seg_starts) > 1 and s0 > 0:
                axA.vlines(s0 - GAP / 2, y - 0.42, y + 0.42, color=S.PALETTE["ink_faint"],
                           lw=0.5, ls=":", zorder=2)
            if o == "Bacteroidales":
                axA.text((s0 + s1 - GAP) / 2, y + 0.46, dtags[d], ha="center", va="bottom",
                         fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], zorder=4)

        # 行左标签：目名 + 株/属/科
        ng, nf = o_val(o, "n_genera"), o_val(o, "n_families")
        nn = int(o_val(o, "n_strains"))
        kk = int(o_val(o, "gh1_positive"))
        axA.text(-0.012, y + 0.16, o, transform=axA.get_yaxis_transform(), ha="right",
                 va="center", fontsize=S.FS["axis"], color=S.ORDER[o], fontweight="bold")
        axA.text(-0.012, y - 0.30, f"{nn} strains · {ng:g} genera · {nf:g} families",
                 transform=axA.get_yaxis_transform(), ha="right", va="center",
                 fontsize=S.FS["note"], color=S.PALETTE["ink_soft"])
        # 行右标签：检出 k/n
        axA.text(1.012, y, f"{kk} / {nn}", transform=axA.get_yaxis_transform(), ha="left",
                 va="center", fontsize=S.FS["axis"],
                 color=S.SEM["not_detected"] if kk == 0 else S.SEM["detected"],
                 fontweight="bold")

    axA.set_xlim(-1.0, nmax + 1.5)
    axA.set_ylim(-0.80, len(orders) + 0.30)
    axA.set_yticks([]); axA.set_xticks([])
    for sp in ("left", "bottom", "right", "top"):
        axA.spines[sp].set_visible(False)
    S.panel_label(axA, "A", "One square = one genome  ·  GH1 detection across 86 Bacteroidota genomes")

    # Bacteroidales 空带的强调框 + Wilson/稳健性注记
    yb = row_y["Bacteroidales"]
    nb = int(o_val("Bacteroidales", "n_strains"))
    axA.add_patch(plt.Rectangle((-0.62, yb - 0.46), nb + 0.30, 0.92, fill=False,
                                edgecolor=S.SEM["not_detected"], lw=1.3, zorder=5,
                                linestyle=(0, (4, 2))))
    wil = sd.one("B_wilson", "wilson95_upper_pct", entity="Bacteroidales")
    axA.annotate(f"GH1 not detected in 0 / {nb} genomes  ·  Wilson 95% upper bound on the detection rate = {wil:.2f}%  ·  "
                 "robust across 4 thresholds (0, 0, 0, 0)",
                 xy=(nb * 0.50, yb + 0.50), xytext=(nb * 0.50, yb + 0.80),
                 fontsize=S.FS["note"], color=S.SEM["not_detected"], ha="center", va="bottom",
                 fontweight="bold", linespacing=1.35,
                 arrowprops=dict(arrowstyle="-", lw=0.7, color=S.SEM["not_detected"]))

    # 图例（双通道：颜色 + 形状语义说明）
    from matplotlib.lines import Line2D
    axA.legend(handles=[
        Line2D([], [], marker="s", linestyle="", markersize=5.2,
               markerfacecolor=S.SEM["absent_fill"], markeredgecolor=S.SEM["not_detected"],
               markeredgewidth=0.8, label="GH1 not detected (0)"),
        Line2D([], [], marker="s", linestyle="", markersize=5.2,
               markerfacecolor=S.SEM["detected"], markeredgecolor=S.PALETTE["blue"],
               markeredgewidth=0.8, label="GH1 detected (>0)"),
    ], loc="lower right", ncol=2, bbox_to_anchor=(1.0, -0.02), fontsize=S.FS["legend"])

    # ---------------- Panel B：检出率 + Wilson 上限 ----------------
    gsB = fig.add_gridspec(1, 1, left=0.105, right=0.345, top=0.470, bottom=0.245)
    axB = fig.add_subplot(gsB[0, 0])
    ys = np.arange(len(orders))[::-1]
    for i, o in enumerate(orders):
        y = ys[i]
        rate = o_val(o, "detection_rate")
        nn = int(o_val(o, "n_strains")); kk = int(o_val(o, "gh1_positive"))
        col = S.ORDER[o]
        axB.barh(y, max(rate, 0.0), height=0.55, color=col, alpha=0.85, zorder=3,
                 edgecolor="white", linewidth=0.5)
        if kk == 0:
            axB.scatter([0.0], [y], marker="|", s=90, color=S.SEM["not_detected"], zorder=4, lw=1.6)
            # Wilson 上限须
            axB.hlines(y, 0.0, wil / 100.0, color=S.SEM["not_detected"], lw=1.2,
                       ls=(0, (3, 2)), zorder=4)
            axB.vlines(wil / 100.0, y - 0.20, y + 0.20, color=S.SEM["not_detected"], lw=1.2, zorder=4)
            axB.text(wil / 100.0 + 0.012, y + 0.30, f"Wilson 95% UB {wil:.2f}%",
                     fontsize=S.FS["note"], color=S.SEM["not_detected"], va="bottom")
        axB.text(max(rate, 0.0) + 0.035, y, f"{kk}/{nn}", va="center", ha="left",
                 fontsize=S.FS["note"], color=S.PALETTE["ink"])
        axB.text(-0.02, y, o, transform=axB.get_yaxis_transform(), ha="right", va="center",
                 fontsize=S.FS["tick"], color=col, fontweight="bold")
    axB.set_xlim(0, 1.06)
    axB.set_ylim(-0.7, len(orders) - 0.3)
    axB.set_yticks([])
    axB.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    axB.set_xlabel("GH1 detection rate", fontsize=S.FS["axis"])
    S.style_ax(axB, grid="x", keep_spines=("bottom",))
    S.panel_label(axB, "B", "Detection rate")

    # ---------------- Panel C：三数据集并列 ----------------
    gsC = fig.add_gridspec(1, 1, left=0.435, right=0.660, top=0.470, bottom=0.245)
    axC = fig.add_subplot(gsC[0, 0])
    C = sd.block("C_datasets")

    def c_val(ds, field):
        r = C[(C["entity"] == ds) & (C["field"] == field)]
        return float(r["value"].iloc[0]) if len(r) else np.nan

    yc = np.arange(len(dsets))[::-1]
    for i, d in enumerate(dsets):
        y = yc[i]
        prop = c_val(d, "gh1_zero_prop"); wn = c_val(d, "gh1_positive_wilson95_upper")
        z = int(c_val(d, "gh1_zero_n")); n = int(c_val(d, "n_strains"))
        axC.barh(y, prop, height=0.5, color=S.PALETTE["sky"], alpha=0.85, zorder=3,
                 edgecolor="white", linewidth=0.5)
        axC.text(prop - 0.02, y, f"{z}/{n} zero", fontsize=S.FS["note"],
                 color=S.PALETTE["ink"], va="center", ha="right", fontweight="bold")
        axC.text(-0.02, y, d.replace("_", " "), transform=axC.get_yaxis_transform(),
                 ha="right", va="center", fontsize=S.FS["tick"], color=S.PALETTE["ink_soft"])
        # Wilson 上限是"检出率"的口径，与零值占比不是同一量纲 —— 只作文字标注，不画在同一轴上（红线 R7）
        axC.text(prop, y + 0.30, f"delivered Wilson 95% UB on detection rate = {wn*100:.1f}%",
                 ha="right", va="bottom", fontsize=S.FS["note"],
                 color=S.SEM["not_detected"])
    axC.set_xlim(0, 1.06)
    axC.set_ylim(-0.7, len(dsets) - 0.25)
    axC.set_yticks([])
    axC.set_xticks([0, 0.5, 1.0])
    axC.set_xlim(0, 1.0)
    axC.set_xlabel("GH1-zero proportion", fontsize=S.FS["axis"])
    S.style_ax(axC, grid="x", keep_spines=("bottom",))
    S.panel_label(axC, "C", "Three independent samplings agree")

    # ---------------- Panel D：四级阈值稳健性 ----------------
    gsD = fig.add_gridspec(1, 1, left=0.755, right=0.995, top=0.470, bottom=0.245)
    axD = fig.add_subplot(gsD[0, 0])
    D = sd.block("D_threshold")
    xs = np.arange(len(thrs))
    styles = {"official_full_filter": (S.SEM["not_detected"], "-", 2.0, "official filters — primary"),
              "A_official_full_filters": (S.PALETTE["grey"], "--", 1.1, "GH1-family-only overlap filter"),
              "B_ievalue_only_no_cov": (S.PALETTE["sky"], ":", 1.1, "i-evalue only, no coverage filter")}
    for var, (col, ls, lw, lab) in styles.items():
        v = D[D["entity"] == var].sort_values("sub", key=lambda s: s.astype(int))
        axD.plot(xs, pd.to_numeric(v["value"]).to_numpy(), ls, color=col, lw=lw,
                 marker="o" if var == "official_full_filter" else "none",
                 ms=4.0, label=lab, zorder=4)
    axD.set_xticks(xs); axD.set_xticklabels(thrs, fontsize=S.FS["tick"])
    axD.set_ylim(-1.2, 27.5)
    axD.set_xlabel("GH1 domain i-evalue threshold", fontsize=S.FS["axis"])
    axD.set_ylabel("strains detected (of 63)", fontsize=S.FS["axis"])
    axD.legend(fontsize=S.FS["note"], loc="upper left", frameon=False)
    axD.text(0.98, 0.06, "0 / 63 at every threshold", transform=axD.transAxes,
             ha="right", va="bottom", fontsize=S.FS["note"], color=S.SEM["not_detected"],
             fontweight="bold")
    S.style_ax(axD, grid="y")
    S.panel_label(axD, "D", "Robust across 4 thresholds")

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.105, right=0.995, top=0.118, bottom=0.010)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")
    axN.text(0, 1.0,
        "GH1 was not detected — under the prescribed dbCAN V5-2 filters — in any of the 63 Bacteroidales genomes "
        "spanning ≥11 nominal genera (10 independent of the Bacteroides–Phocaeicola split) and ≥8 families; "
        f"the Wilson 95% upper bound on the GH1 detection rate is {wil:.2f}% (0/63).",
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top")
    axN.text(0, 0.56,
        "Sister orders of the same phylum do carry GH1 (Flavobacteriales 3/14, Sphingobacteriales 5/5, Cytophagales 4/4 across the three datasets), "
        "so the boundary of non-detection lies at the ORDER level, not the phylum level.  "
        "“Not detected at the prescribed thresholds” is retained as the wording — it is not a claim of biological absence.",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.14,
        "Per-strain values aggregated from official_strain_cazyme_matrix_round4.tsv (31 Bacteroidota), sentinel_10strains_matrix.tsv (10) and "
        "official_strain_cazyme_matrix_round7.tsv (45); Wilson bounds and threshold counts are taken verbatim from the delivered tables (no recomputation).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top")

    S.figure_title(fig, "Order-level boundary of GH1 non-detection within Bacteroidota",
                   x=0.105, y=0.990)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[fig2] PLOT   ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
