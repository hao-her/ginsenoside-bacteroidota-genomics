#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_suppfig2.py — Suppl Fig 2 系统发育树 + GH1 映射（环形线稿，1200 dpi）
=====================================================================
【任务书要求】树 + 旁侧热图；必须标注 "tree covers 114 strains; 55 newly sampled
strains not included (no sequences)"（该整图注在返工轮移至 Word 图注，图面不再嵌入）。

【数据来源说明】任务书指的 outputs_round5/trees/bacteria_annotated_*.nwk 未随附；
本图使用随附的 outputs_round8/stats/bacteria_midpoint.nwk（114 个 tip，经核对与
official_strain_cazyme_matrix_round4.tsv 的 114 株细菌 accession 完全一致），
分类与 GH 计数由 round4 矩阵按 accession 关联（逐行 provenance）。

【美化轮设计思路 —— 环形（fan）布局】
旧版为纵向矩形树 + 旁侧热图，纵高 212 mm、字体被压小、顶部 colorbar 刻度与旋转
的家族列标签重叠（见返工轮报告 §5.1）。本轮改为环形布局：
  · 中心 → 外：branch-length 环形谱系树（midpoint 根），tip 以门色小点标记；
  · 同心环带（内→外）：门色环 → GH1 二值环（朱红空 = 未检出 / 蓝实 = 检出）
    → 10 家族拷贝数热图环（单色相明度渐变，天然色盲安全）；
  · 外圈：属级连续片段弧括号 + 大字号属名（树序下恰好 8 个连续片段）；
  · 顶部缺口：12 条环带的引出式标注（callout），彻底避免旧版列标签重叠；
  · 中心留白：门 / GH1 图例 + colorbar；
  · Bacteroidota 31 株整段：GH1 环外加朱红强调弧 + 注记（核心信息视觉锚点）。
环形布局在同一画幅内让每个标签的可用弧长更大 → 字号更大、占版更小。
所有绘制数值与旧版逐一相同（同一 source_data.tsv + 同一 nwk），零重算、零改数。

运行：python plot_suppfig2.py / --plot-only
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Wedge
from matplotlib.lines import Line2D

import shared_style as S

FIG_NAME = "suppfig2_phylogeny_gh1_mapping"
OUT_DIR = S.SUPP_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = 190.0, 190.0        # 方形画幅（环形）
KIND = "line"

FAMS = ["GH1", "GH2", "GH3", "GH27", "GH36", "GH51", "GH54", "GH78", "GH79", "GH106"]


# ======================================================================
# newick 解析（极简递归下降）
# ======================================================================
def parse_newick(text: str):
    s = text.strip()
    pos = 0

    def node():
        nonlocal pos
        nd = {"ch": [], "name": None, "bl": 0.0}
        if s[pos] == "(":
            pos += 1
            while True:
                nd["ch"].append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^():,;]+", s[pos:])
        if m:
            nd["name"] = m.group(0).strip()
            pos += m.end()
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m2 = re.match(r"[0-9.eE+-]+", s[pos:])
            nd["bl"] = float(m2.group(0))
            pos += m2.end()
        return nd

    root = node()
    return root


def tip_order(root):
    order = []

    def rec(n):
        if not n["ch"]:
            order.append(n["name"])
            return
        for c in n["ch"]:
            rec(c)
    rec(root)
    return order


# ======================================================================
# BUILD（与上一轮逐字一致；plot-only 模式不调用，不触碰上游计算表）
# ======================================================================
def build() -> S.SourceTable:
    st = S.SourceTable("Suppl Fig 2 — bacterial phylogeny (114 strains) with GH-family mapping")
    P_tree, P_m4 = S.rel(S.SRC["tree_nwk"]), S.rel(S.SRC["matrix_r4"])

    root = parse_newick(Path(S.SRC["tree_nwk"]).read_text())
    order = tip_order(root)
    m4 = S.read_upstream("matrix_r4").set_index("accession")

    for i, acc in enumerate(order):
        r = m4.loc[acc]
        st.add("tips", acc, str(i), "phylum", str(r["phylum"]), prov=P_m4)
        st.add("tips", acc, str(i), "genus", str(r["genus"]), prov=P_m4)
        st.add("tips", acc, str(i), "gh1", str(int(r["GH1"])), prov=P_m4)
        for f in FAMS:
            st.add("heatmap", acc, str(i), f, str(int(r[f])), prov=P_m4)
    st.add("meta", "n_tips", "", "value", str(len(order)), prov=P_tree)
    st.add("meta", "n_not_in_tree", "", "value", "55",
           text="round7 45 + sentinel 10 newly sampled strains — no sequences available", prov="handoff §4 / round8 B5")
    d = S.read_upstream("phylo_D"); lam = S.read_upstream("phylo_lambda")
    st.add("meta", "D", "", "value", str(d["D"].iloc[0]), prov=S.rel(S.SRC["phylo_D"]))
    for _, r in lam.iterrows():
        st.add("meta", f"lambda_{r['family']}", "", "value", str(r["lambda_ML"]), prov=S.rel(S.SRC["phylo_lambda"]))
    for i, f in enumerate(FAMS):
        st.add("meta", "family_order", str(i), "name", text=f, prov="10-family panel")

    st.write(SOURCE_DATA)
    print(f"[supp2] BUILD -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# 环形几何小工具
# ======================================================================
GAP_DEG = 48.0                       # 顶部缺口（放环带引出标注）
THETA_TOP = 90.0                     # 缺口中心角

def _polar(r, theta_deg):
    a = np.radians(theta_deg)
    return r * np.cos(a), r * np.sin(a)


# ======================================================================
# PLOT（环形重绘）
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    T = sd.block("tips")
    H = sd.block("heatmap")
    n = int(sd.one("meta", "value", entity="n_tips"))
    fams = [t for t in (sd.one_text("meta", "name", entity="family_order", sub=str(i))
                        for i in range(len(FAMS))) if t]

    def _series(field, numeric=False):
        v = T[T["field"] == field].sort_values("sub", key=lambda s: s.astype(int))["value"]
        return pd.to_numeric(v).tolist() if numeric else v.tolist()
    accs = T[T["field"] == "phylum"].sort_values("sub", key=lambda s: s.astype(int))["entity"].tolist()
    phylum = _series("phylum")
    genus = _series("genus")
    gh1 = _series("gh1", numeric=True)

    Hm = (H.pivot_table(index="sub", columns="field", values="value", aggfunc="first")
          .reindex([str(i) for i in range(n)])[fams].apply(pd.to_numeric))

    root = parse_newick(Path(S.SRC["tree_nwk"]).read_text())
    yidx = {a: i for i, a in enumerate(accs)}

    # ---- 树布局：x = 累积枝长，y = tip 行号（内部节点取子均值）----
    coords = {}

    def layout(nd, x0):
        x = x0 + nd["bl"]
        if not nd["ch"]:
            y = yidx[nd["name"]]
            coords[id(nd)] = (x, y)
            return x, y
        ys = []
        for c in nd["ch"]:
            _, cy = layout(c, x)
            ys.append(cy)
        y = float(np.mean(ys))
        coords[id(nd)] = (x, y)
        return x, y

    layout(root, 0.0)
    xmax = max(v[0] for v in coords.values()) or 1.0

    # ---- 半径 / 角度映射 ----
    SWEEP = 360.0 - GAP_DEG
    dt = SWEEP / n                                  # 每 tip 角宽（度）
    a_start = THETA_TOP - GAP_DEG / 2.0             # 缺口右缘，顺时针铺开

    def theta_of(y):                                # tip 行号(可为浮点) → 角度
        return a_start - (y + 0.5) * dt

    R_INNER, R_TREE = 0.34, 0.605                   # 树环
    def r_of(x):
        return R_INNER + (x / xmax) * (R_TREE - R_INNER)

    # 同心环带半径（内缘, 外缘）
    R_PH = (0.618, 0.643)                           # 门色环
    R_G1 = (0.650, 0.691)                           # GH1 二值环（加粗 = 主角）
    HM0, PITCH, HMW = 0.703, 0.0272, 0.0232         # 热图环起点/步距/环宽
    R_HM = [(HM0 + k * PITCH, HM0 + k * PITCH + HMW) for k in range(len(fams))]
    R_BR = R_HM[-1][1] + 0.010                      # 属括号弧半径
    R_GLAB = R_BR + 0.012                           # 属名起始半径

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_aspect("equal")
    ax.set_xlim(-1.62, 1.62)
    ax.set_ylim(-1.62, 1.62)
    ax.axis("off")

    ink_soft, ink_faint = S.PALETTE["ink_soft"], S.PALETTE["ink_faint"]

    # ---- 环形树 ----
    def draw(nd):
        rx, ry = r_of(coords[id(nd)][0]), coords[id(nd)][1]
        if nd["ch"]:
            cths = [theta_of(coords[id(c)][1]) for c in nd["ch"]]
            aa = np.linspace(min(cths), max(cths), 60)
            xs, ys = _polar(rx, aa)
            ax.plot(xs, ys, color=ink_soft, lw=0.5, solid_capstyle="round", zorder=2)
            for c in nd["ch"]:
                cth = theta_of(coords[id(c)][1])
                cx = r_of(coords[id(c)][0])
                x0, y0 = _polar(rx, cth); x1, y1 = _polar(cx, cth)
                ax.plot([x0, x1], [y0, y1], color=ink_soft, lw=0.5,
                        solid_capstyle="round", zorder=2)
            for c in nd["ch"]:
                draw(c)
        else:
            th = theta_of(ry)
            # tip → 环带的浅色对齐引导线
            gx0, gy0 = _polar(r_of(coords[id(nd)][0]), th)
            gx1, gy1 = _polar(R_PH[0] - 0.004, th)
            ax.plot([gx0, gx1], [gy0, gy1], color=S.PALETTE["grid"], lw=0.30, zorder=1)
            col = S.PHYLUM.get(phylum[int(ry)], S.PALETTE["grey"])
            tx, ty = _polar(r_of(coords[id(nd)][0]), th)
            ax.scatter([tx], [ty], marker="o", s=4.5, color=col, edgecolors="white",
                       linewidths=0.20, zorder=4)
    draw(root)

    # ---- 环带绘制器 ----
    def ring(i, r0, r1, facecolor, edgecolor="none", lw=0.0, zorder=3):
        th = theta_of(i)
        Wedge_ = Wedge((0, 0), r1, th - dt / 2.0, th + dt / 2.0, width=r1 - r0,
                       facecolor=facecolor, edgecolor=edgecolor, lw=lw, zorder=zorder)
        ax.add_patch(Wedge_)

    # 数据扇区角度范围（用于画环边界弧）
    a_hi = theta_of(0) + dt / 2.0
    a_lo = theta_of(n - 1) - dt / 2.0

    def bound_arc(r, color=S.PALETTE["axis"], lw=0.18, zorder=2):
        aa = np.linspace(a_lo, a_hi, 400)
        xs, ys = _polar(r, aa)
        ax.plot(xs, ys, color=color, lw=lw, zorder=zorder, solid_capstyle="butt")

    # 门色环
    for i, ph in enumerate(phylum):
        ring(i, R_PH[0], R_PH[1], S.PHYLUM.get(ph, S.PALETTE["grey"]))
    # GH1 二值环
    for i, g in enumerate(gh1):
        if g > 0:
            ring(i, R_G1[0], R_G1[1], S.SEM["detected"])
        else:
            ring(i, R_G1[0], R_G1[1], S.SEM["absent_fill"],
                 edgecolor=S.SEM["not_detected"], lw=0.28)
    # 10 家族热图环
    cmap = LinearSegmentedColormap.from_list("seq_oi", S.SEQ_RAMP)
    vmax = float(Hm.to_numpy().max())
    for k, f in enumerate(fams):
        vals = Hm[f].to_numpy(dtype=float)
        r0, r1 = R_HM[k]
        for i in range(n):
            ring(i, r0, r1, cmap(vals[i] / vmax))
    # 淡色环边界弧（让空值环带仍具结构感）
    bound_arc(R_PH[0]); bound_arc(R_PH[1])
    bound_arc(R_G1[0], color=S.SEM["not_detected"], lw=0.22)
    bound_arc(R_G1[1], color=S.SEM["not_detected"], lw=0.22)
    for k in range(len(fams)):
        bound_arc(R_HM[k][1], color=S.PALETTE["grid"], lw=0.16)
    bound_arc(R_HM[0][0], color=S.PALETTE["grid"], lw=0.16)

    # ---- Bacteroidota 整段：GH1 环外朱红强调弧 + 注记 ----
    bact = [i for i, ph in enumerate(phylum) if ph == "Bacteroidota"]
    if bact:
        b0, b1 = min(bact), max(bact)
        t_lo = theta_of(b1) - dt / 2.0
        t_hi = theta_of(b0) + dt / 2.0
        aa = np.linspace(t_lo, t_hi, 200)
        rr = R_G1[1] + 0.012
        xs, ys = _polar(rr, aa)
        ax.plot(xs, ys, color=S.SEM["not_detected"], lw=1.8, solid_capstyle="round", zorder=6)
        # 两端径向小帽
        for tt in (t_lo, t_hi):
            x0, y0 = _polar(R_G1[0] - 0.004, tt); x1, y1 = _polar(rr, tt)
            ax.plot([x0, x1], [y0, y1], color=S.SEM["not_detected"], lw=1.3, zorder=6)
        # 注记：置于左上角留白，细引线指向强调弧靠 Parabacteroides 端
        ax_pt = _polar(rr, t_lo + 8)
        txt_x, txt_y = -1.42, 1.20
        ax.annotate("Bacteroidota — 31 strains\nGH1 = 0 in every genome",
                    xy=ax_pt, xytext=(txt_x, txt_y),
                    fontsize=S.FS["annot"], color=S.SEM["not_detected"], fontweight="bold",
                    ha="left", va="center", linespacing=1.35, zorder=9,
                    arrowprops=dict(arrowstyle="-", lw=0.7, color=S.SEM["not_detected"],
                                    connectionstyle="arc3,rad=-0.15", shrinkA=2, shrinkB=3))

    # ---- 属级连续片段弧括号 + 属名 ----
    runs = []
    for i, g in enumerate(genus):
        if not runs or runs[-1][0] != g:
            runs.append([g, i, i])
        else:
            runs[-1][2] = i
    # 属 → 门（用于标签着色）
    g2p = {}
    for i, g in enumerate(genus):
        g2p.setdefault(g, phylum[i])
    for g, a, b in runs:
        t_lo = theta_of(b) - dt / 2.0 + 0.6
        t_hi = theta_of(a) + dt / 2.0 - 0.6
        aa = np.linspace(t_lo, t_hi, 80)
        xs, ys = _polar(R_BR, aa)
        col = S.PHYLUM.get(g2p[g], ink_soft)
        ax.plot(xs, ys, color=col, lw=1.4, solid_capstyle="round", zorder=5)
        # 括号两端小径向帽
        for tt in (t_lo, t_hi):
            x0, y0 = _polar(R_BR, tt); x1, y1 = _polar(R_BR - 0.008, tt)
            ax.plot([x0, x1], [y0, y1], color=col, lw=1.4, zorder=5)
        # 属名（径向朝外）
        tmid = 0.5 * (t_lo + t_hi)
        lx, ly = _polar(R_GLAB, tmid)
        rot = tmid
        ha = "left"
        if 90 < (tmid % 360) < 270:
            rot += 180
            ha = "right"
        ax.text(lx, ly, f"{g}  ({b - a + 1})", fontsize=S.FS["axis"] + 1.0,
                color=col, fontweight="bold", rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", zorder=8)

    # ---- 顶部缺口：12 条环带的引出式标注（callout ladder）----
    ring_specs = [("phylum", 0.5 * (R_PH[0] + R_PH[1]), ink_soft, "normal")]
    ring_specs.append(("GH1 detected", 0.5 * (R_G1[0] + R_G1[1]), S.SEM["not_detected"], "bold"))
    for k, f in enumerate(fams):
        ring_specs.append((f, 0.5 * (R_HM[k][0] + R_HM[k][1]),
                           S.FAMILY.get(f, ink_soft),
                           "bold" if f in ("GH1", "GH2", "GH3") else "normal"))
    m = len(ring_specs)
    # 在缺口角度内均匀分配标注角，径向朝外，末端带引线到真实环半径
    lab_lo = THETA_TOP - GAP_DEG / 2.0 + 3.5
    lab_hi = THETA_TOP + GAP_DEG / 2.0 - 3.5
    lab_angles = np.linspace(lab_hi, lab_lo, m)      # 由左到右
    R_TICK = R_HM[-1][1] + 0.028
    R_TXT = R_TICK + 0.016
    for (name, r_mid, col, weight), ang in zip(ring_specs, lab_angles):
        # 引线：真实环半径(在缺口右缘) → 标注角的刻度点
        x_ring, y_ring = _polar(r_mid, THETA_TOP - GAP_DEG / 2.0 + 0.5)
        x_tick, y_tick = _polar(R_TICK, ang)
        ax.plot([x_ring, x_tick], [y_ring, y_tick], color=col, lw=0.5, alpha=0.8, zorder=6)
        ax.scatter([x_ring], [y_ring], s=5.0, color=col, zorder=7,
                   edgecolors="white", linewidths=0.2)
        lx, ly = _polar(R_TXT, ang)
        rot = ang
        ha = "left"
        if 90 < (ang % 360) < 270:
            rot += 180; ha = "right"
        ax.text(lx, ly, name, fontsize=S.FS["tick"], color=col, fontweight=weight,
                rotation=rot, rotation_mode="anchor", ha=ha, va="center", zorder=8)

    # ---- 中心留白：门 / GH1 图例 + colorbar ----
    # 门图例
    ph_present = [p for p in S.PHYLUM_ORDER if p in set(phylum)]
    handles = [Line2D([], [], marker="o", linestyle="", markersize=5.2,
                      markerfacecolor=S.PHYLUM[p], markeredgecolor="white",
                      markeredgewidth=0.3, label=p) for p in ph_present]
    handles += [
        Line2D([], [], marker="s", linestyle="", markersize=5.2,
               markerfacecolor=S.SEM["detected"], markeredgecolor="none", label="GH1 detected (>0)"),
        Line2D([], [], marker="s", linestyle="", markersize=5.2,
               markerfacecolor=S.SEM["absent_fill"], markeredgecolor=S.SEM["not_detected"],
               markeredgewidth=0.7, label="GH1 not detected (0)"),
    ]
    leg = ax.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.55),
                    bbox_transform=fig.transFigure, frameon=False,
                    fontsize=S.FS["legend"], labelspacing=0.5, handletextpad=0.5,
                    borderaxespad=0.0)
    for txt in leg.get_texts():
        txt.set_color(S.PALETTE["ink"])

    # colorbar（水平，置于中心图例下方）
    cbar_ax = fig.add_axes([0.440, 0.438, 0.12, 0.0095])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax)),
                      cax=cbar_ax, orientation="horizontal")
    cb.set_ticks([0, vmax / 2, vmax])
    cb.set_ticklabels(["0", f"{vmax/2:g}", f"{vmax:g}"])
    cb.ax.tick_params(labelsize=S.FS["note"], length=1.5, pad=1.0)
    cb.outline.set_visible(False)
    fig.text(0.5, 0.431, "GH copy number / genome", fontsize=S.FS["note"],
             color=ink_soft, ha="center", va="top")

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[supp2] PLOT  ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
