#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
shared_style.py — 全局样式参数（绘图轮唯一样式真源）
=====================================================
所有 plot_*.py 都从这里取：字体、字号、配色、线宽、图幅尺寸、dpi、
坐标轴装饰、三格式保存器、以及 source_data.tsv 的读写器。

修改本文件即可全局生效（绘图指令「脚本交付要求」第 4 条）。

配色原则
--------
全部色值取自 Okabe–Ito (2008) 色盲友好色板，禁用红绿对比。
语义在全套 8 张图中保持一致（见 SEM / PHYLUM / FAMILY / ORDER 字典），
因此跨图不会出现"同一颜色代表不同含义"的割裂感。
色盲验证见 palette/ 与 scripts/make_palette_sheet.py。

数据溯源原则
------------
脚本内不写死任何绘图数值。每个 plot_*.py 分两个阶段：
  1) BUILD  —— 从上游项目表读取，装配本图的 <fig>_source_data.tsv
  2) PLOT   —— 只从 <fig>_source_data.tsv 读数并绘图
SourceTable / read_source 就是这两阶段之间的唯一接口。
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # 无交互后端（绘图指令第 2 条：不弹窗）

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

# ======================================================================
# 0. 路径（全部集中在顶部变量，函数内不硬编码 —— 绘图指令第 2 条）
# ======================================================================
HERE       = Path(__file__).resolve().parent          # .../outputs_figures/scripts
FIGROOT    = HERE.parent                              # .../outputs_figures
MAIN_DIR   = FIGROOT / "main"
SUPP_DIR   = FIGROOT / "supplementary"
PALETTE_DIR = FIGROOT / "palette"

# 上游只读输入目录（round4/5/7/8/9 交付件）。
# 默认按本交付包相对位置推导（scripts/ 的上三级 sibling = input/），不含硬编码绝对路径；
# 可用环境变量 FIG_INPUT 覆盖。
INPUT = Path(os.environ.get("FIG_INPUT", str(HERE.parents[2] / "input")))

# 上游文件（相对 INPUT）。集中登记 = 单一事实来源，便于审计。
SRC = {
    # --- Fig 1 ---
    "matrix_r4":        INPUT / "official_strain_cazyme_matrix_round4.tsv",
    "phylum_dist":      INPUT / "phylum_family_distribution.tsv",
    "genus_dist":       INPUT / "genus_family_distribution.tsv",
    "genus_kw":         INPUT / "outputs_round8/tables/t3a_genus_kw.tsv",
    "strain_kw":        INPUT / "kruskal_epsilon2.tsv",
    "phylo_D":          INPUT / "outputs_round8/tables/t3b_D_statistic.tsv",
    "phylo_lambda":     INPUT / "t3b_continuous_signal.tsv",
    "genus_medians":    INPUT / "t3a_genus_medians.tsv",
    # --- Fig 2 ---
    "matrix_r7":        INPUT / "official_strain_cazyme_matrix_round7.tsv",
    "matrix_sentinel":  INPUT / "sentinel_10strains_matrix.tsv",
    "three_datasets":   INPUT / "task3_2_three_datasets_parallel.tsv",
    "t4_prov":          INPUT / "tables/t4_number_provenance.tsv",
    "gh1_curve":        INPUT / "t5_gh1_threshold_curve.tsv",
    "gh23_curve":       INPUT / "t5_gh23_control_curves.tsv",
    # --- Fig 3 ---
    "condnull":         INPUT / "tables/t2_condnull_aggregate.tsv",
    "uniformnull":      INPUT / "outputs_round8/tables/t1_nullmodel_aggregate.tsv",
    "unbiased45":       INPUT / "tables/t3_null45_aggregate.tsv",
    "unbiased45_layer": INPUT / "tables/t3_null45_perlayer.tsv",
    "unbiased31":       INPUT / "tables/t3_null31all_aggregate.tsv",
    "unbiased31_genus": INPUT / "tables/t3_null31all_pergenus.tsv",
    "fisher_sidechain": INPUT / "outputs_round8/tables/t2c_fisher_tests.tsv",
    "fisher_loo":       INPUT / "outputs_round8/tables/t2c_loo_sensitivity.tsv",
    "occupancy":        INPUT / "outputs_round8/tables/t2a_occupancy_aggregate.tsv",
    # --- Fig 4 ---
    # ⚠ 任务书指的 学生 TableA_cazy_activities.xlsx 未随附（见 report_figures.md 数据缺口 D1）。
    #   Activities/Characterized 数取自随附的《绘图轮全景任务表.MD》§2 表（同一项目的已定稿数值）。
    "handoff_md":       INPUT / "绘图轮全景任务表.MD",
    "power_curve":      INPUT / "outputs_round8/figures/fig3_source_data.tsv",
    "power_grid":       INPUT / "outputs_round8/tables/t6_scenario_grid.tsv",
    "power_mc":         INPUT / "outputs_round8/tables/t6_montecarlo.tsv",
    # --- Fig 5 ---
    "gh3_cand":         INPUT / "tables/t1_gh3_candidates_all.tsv",
    "gh5_cand":         INPUT / "tables/t1_gh5_candidates_all.tsv",
    "gh66_inv":         INPUT / "tables/t1_gh66_inventory.tsv",
    "ranked":           INPUT / "tables/t1_ranked_passers.tsv",
    "shortlist":        INPUT / "tables/t1_wetlab_shortlist_clusters.tsv",
    "pool_summary":     INPUT / "stats/t1_pool_summary.json",
    # --- Suppl Fig 1 ---
    "frag_inventory":   INPUT / "t5_gh1_fragment_inventory.tsv",
    # --- Suppl Fig 2 ---
    "tree_nwk":         INPUT / "outputs_round8/stats/bacteria_midpoint.nwk",
    # --- Suppl Fig 3 ---
    "strain45":         INPUT / "tables/t3_strain_summary_45.tsv",
    "pairs45":          INPUT / "tables/t3_pairs_45.tsv",
}


# ======================================================================
# 1. 配色 —— Okabe–Ito 色盲友好色板（全套图统一）
# ======================================================================
PALETTE = {
    # 墨色与结构色（非数据色）
    "ink":        "#1A1A1A",   # 主文字
    "ink_soft":   "#5C5C5C",   # 次级文字
    "ink_faint":  "#8C8C8C",   # 弱化文字
    "grid":       "#E8E8E8",   # 网格
    "axis":       "#B8B8B8",   # 轴线
    "paper":      "#FFFFFF",
    "panel":      "#F6F6F6",   # 面板底
    "panel_edge": "#DCDCDC",

    # Okabe–Ito 数据色（7 色）
    "blue":       "#0072B2",
    "vermillion": "#D55E00",
    "green":      "#009E73",
    "orange":     "#E69F00",
    "purple":     "#CC79A7",
    "sky":        "#56B4E9",
    "grey":       "#646464",
    "yellow":     "#F0E442",   # Okabe-Ito yellow   # 中性"无信号"色；取较深值以在 deuteranopia 下与 purple/green 保持 ΔE≥15

    # 派生淡色（同色相，仅改明度 → 不新增色盲风险）
    "blue_l":     "#B8D9EA",
    "verm_l":     "#F6CDB6",
    "green_l":    "#B4E0D2",
    "orange_l":   "#F7DFAE",
    "purple_l":   "#EED3E2",
    "sky_l":      "#CDE7F6",
    "grey_l":     "#E0E0E0",
}

# --- 语义映射：一个颜色在全套图中只代表一个含义 -----------------------
SEM = {
    # 检出 / 未检出（Fig 2、Suppl 1、Suppl 2 的核心对立）
    "detected":       PALETTE["blue"],        # GH1 检出（阳性）
    "not_detected":   PALETTE["vermillion"],  # GH1 未检出（醒目强调，任务书指定）
    "absent_fill":    PALETTE["verm_l"],
    "control":        PALETTE["green"],       # 阳性对照 GH2/GH3
    "context":        PALETTE["grey"],        # 仅作背景/上下文，不作结论

    # PUL 双口径（Fig 3、Suppl 3）——两口径必须视觉分离（红线 R7）
    "conditional":    PALETTE["vermillion"],  # 条件化口径 = 上界
    "unbiased":       PALETTE["blue"],        # 无偏口径 = 主估计
    "significant":    PALETTE["blue"],
    "not_sig":        PALETTE["grey"],

    # 家族角色（Fig 4A、Fig 1）
    "role_core":      PALETTE["blue"],        # 核心 GH1/GH3
    "role_secondary": PALETTE["orange"],      # 次级 GH2
    "role_sidechain": PALETTE["purple"],      # 共现侧链糖苷酶
    "role_nodata":    PALETTE["grey_l"],      # 数据缺失（留空，不插补 —— 红线 R1）

    # 候选筛选（Fig 5）
    "pass":           PALETTE["blue"],
    "rescued":        PALETTE["sky"],
    "fail":           PALETTE["grey"],
    "not_computable": PALETTE["grey_l"],
    "highlight":      PALETTE["vermillion"],
}

# 门级配色（Fig 1 / Suppl 2）
PHYLUM = {
    "Bacteroidota":     PALETTE["blue"],
    "Bacillota":        PALETTE["orange"],
    "Ascomycota":       PALETTE["green"],
    "Actinomycetota":   PALETTE["purple"],
}
PHYLUM_ORDER = ["Bacteroidota", "Bacillota", "Actinomycetota", "Ascomycota"]

# GH 家族配色（Fig 1 三面板 / Suppl 2 热图边条）
FAMILY = {
    "GH1":   PALETTE["blue"],        # 核心
    "GH3":   PALETTE["green"],       # 核心
    "GH2":   PALETTE["orange"],      # 次级
    "GH27":  PALETTE["purple"],
    "GH36":  PALETTE["sky"],
    "GH51":  PALETTE["vermillion"],
    "GH54":  PALETTE["grey"],
    "GH78":  PALETTE["ink_soft"],
    "GH79":  PALETTE["orange_l"],
    "GH106": PALETTE["purple_l"],
}
FAMILY_ROLE = {
    "GH1": "core", "GH3": "core", "GH2": "secondary",
    "GH27": "sidechain", "GH36": "sidechain", "GH51": "sidechain",
    "GH78": "sidechain", "GH106": "sidechain",
    "GH54": "unresolved", "GH79": "unresolved",
}

# 目级配色（Fig 2 / Suppl 3）
ORDER = {
    "Bacteroidales":      PALETTE["vermillion"],   # 缺失带 —— 视觉锚点
    "Flavobacteriales":   PALETTE["blue"],
    "Sphingobacteriales": PALETTE["green"],
    "Cytophagales":       PALETTE["orange"],
}
ORDER_ROW = ["Bacteroidales", "Flavobacteriales", "Sphingobacteriales", "Cytophagales"]

# 属级配色（Fig 5 散点；与位置编码冗余 → 即使两色混淆也不丢信息）
GENUS = {
    "Bacteroides":      PALETTE["blue"],
    "Parabacteroides":  PALETTE["vermillion"],
    "Phocaeicola":      PALETTE["sky"],
    "Barnesiella":      PALETTE["green"],
    "Muribaculum":      PALETTE["orange"],
    "Prevotella":       PALETTE["purple"],
    # n ≤ 4 的稀有属共用中性灰：其身份由行标签文字承担（颜色仅冗余通道），
    # 避免为 11 个属硬凑 11 个色盲可分色而牺牲可分辨性。
    "Alistipes":        PALETTE["yellow"],
    "Tannerella":       PALETTE["grey"],
    "Odoribacter":      PALETTE["grey"],
    "Rikenella":        PALETTE["grey"],
    "Porphyromonas":    PALETTE["grey"],
}

# 顺序色带（Suppl 2 热图；单色相明度递进 → 天然色盲安全）
SEQ_RAMP = ["#FFFFFF", "#E4EEF5", "#B8D9EA", "#7FB6D8", "#3D87BF", "#0072B2", "#004E7C"]


# ======================================================================
# 2. 字体 / 字号 / 线宽 / 图幅 / dpi
# ======================================================================
FONT_STACK = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
# Liberation Sans 与 Arial 度量兼容；本机无 Arial 时自动落到它，字形宽度一致。

# 字号（pt）——红线：正文 ≥7pt，坐标轴标签 ≥8pt
FS = {
    "fig_title":  10.0,
    "panel":       9.0,   # 面板字母 + 标题
    "axis":        8.0,   # 坐标轴标签（≥8pt 硬要求）
    "tick":        7.5,
    "annot":       7.5,
    "note":        7.0,   # 正文最小值（≥7pt 硬要求）
    "legend":      7.5,
    "stat":        7.5,
}

LW = {
    "axis":    0.7,
    "grid":    0.4,
    "box":     0.8,
    "curve":   1.4,
    "curve_thin": 0.9,
    "accent":  1.8,
    "halo":    0.0,
}

# 图幅（mm）——期刊三档：单栏 88 / 1.5 栏 120 / 双栏 183
MM = {"single": 88.0, "mid": 120.0, "double": 183.0}
DPI_MAIN = 600     # 主图
DPI_LINE = 1200    # 线稿/文字图

# 期刊图常见的默认边距比例
PAD = 0.02


def mm(v: float) -> float:
    """mm → inch"""
    return v / 25.4


def figsize(w_mm: float, h_mm: float) -> tuple[float, float]:
    return (mm(w_mm), mm(h_mm))


# ======================================================================
# 3. rcParams 全局装配
# ======================================================================
def apply_rc() -> None:
    """一次性设定全局 rcParams。所有脚本入口调用。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    chosen = next((f for f in FONT_STACK if f in available), FONT_STACK[-1])

    matplotlib.rcParams.update({
        "font.family":      "sans-serif",
        "font.sans-serif":  [chosen] + FONT_STACK,
        "font.size":        FS["note"],

        "axes.titlesize":   FS["panel"],
        "axes.labelsize":   FS["axis"],
        "axes.labelcolor":  PALETTE["ink"],
        "axes.edgecolor":   PALETTE["axis"],
        "axes.linewidth":   LW["axis"],
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.grid":        False,
        "axes.axisbelow":   True,

        "xtick.labelsize":  FS["tick"],
        "ytick.labelsize":  FS["tick"],
        "xtick.color":      PALETTE["ink_soft"],
        "ytick.color":      PALETTE["ink_soft"],
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": LW["axis"],
        "ytick.major.width": LW["axis"],
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.direction":  "out",
        "ytick.direction":  "out",

        "legend.fontsize":  FS["legend"],
        "legend.frameon":   False,
        "legend.labelspacing": 0.35,
        "legend.handlelength": 1.4,
        "legend.borderaxespad": 0.2,

        "lines.linewidth":  LW["curve"],
        "lines.markersize": 3.2,
        "patch.linewidth":  LW["box"],
        "patch.edgecolor":  PALETTE["ink"],
        "patch.force_edgecolor": False,

        "text.color":       PALETTE["ink"],
        "savefig.facecolor": "white",
        "savefig.edgecolor": "none",
        "savefig.dpi":      DPI_MAIN,
        "savefig.bbox":     "tight",
        "savefig.pad_inches": PAD,

        # 矢量输出：SVG 转路径（无字体依赖，任何机器渲染一致）；
        #            PDF/PS 用 TrueType(42) 内嵌，便于 Illustrator 编辑
        "svg.fonttype": "path",
        "pdf.fonttype": 42,
        "ps.fonttype":  42,

        "figure.dpi":   150,
        "figure.constrained_layout.use": False,
    })
    return chosen


# ======================================================================
# 4. 坐标轴 / 面板装饰
# ======================================================================
def style_ax(ax, grid: str = "y", keep_spines=("left", "bottom"),
             grid_color=None, grid_lw=None):
    """统一坐标轴外观：去脊线、加淡网格。grid ∈ {'x','y','both',None}"""
    gc = grid_color or PALETTE["grid"]
    gl = grid_lw or LW["grid"]
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep_spines)
    if grid in ("y", "both"):
        ax.yaxis.grid(True, color=gc, lw=gl, zorder=0)
    if grid in ("x", "both"):
        ax.xaxis.grid(True, color=gc, lw=gl, zorder=0)
    ax.set_axisbelow(True)
    return ax


def panel_label(ax, letter: str, title: str = "", dx: float = 0.0,
                dy: float = 0.0, color=None, size=None, weight="bold"):
    """左上角面板字母（粗体）+ 紧随其后的常规字重标题。"""
    c = color or PALETTE["ink"]
    sz = size or FS["panel"]
    ax.text(dx - 0.005, 1.025 + dy, letter, transform=ax.transAxes,
            fontsize=sz + 1.0, fontweight=weight, color=c,
            va="bottom", ha="left", zorder=10)
    if title:
        ax.text(dx + 0.028, 1.025 + dy, title, transform=ax.transAxes,
                fontsize=sz - 0.5, color=PALETTE["ink_soft"],
                va="bottom", ha="left", zorder=10)
    return ax


def note(ax, text: str, x: float = 0.0, y: float = -0.20, size=None,
         color=None, ha="left", va="top", style="normal", **kw):
    """轴外注记（图注性文字，必须如实呈现的科学信息走这里）。"""
    return ax.text(x, y, text, transform=ax.transAxes,
                   fontsize=size or FS["note"], color=color or PALETTE["ink_soft"],
                   ha=ha, va=va, style=style, linespacing=1.45, **kw)


def figure_note(fig, text: str, x: float = 0.0, y: float = 0.005,
                size=None, color=None, ha="left", va="bottom", **kw):
    """整图底部的强制科学注记（Wilson 上限 / 双口径 / 树覆盖范围等）。"""
    return fig.text(x, y, text, fontsize=size or FS["note"],
                    color=color or PALETTE["ink_soft"], ha=ha, va=va,
                    linespacing=1.5, **kw)


def figure_title(fig, text: str, x: float = 0.0, y: float = 0.995,
                 size=None, color=None, ha="left", va="top", **kw):
    return fig.text(x, y, text, fontsize=size or FS["fig_title"],
                    color=color or PALETTE["ink"], ha=ha, va=va,
                    fontweight="bold", linespacing=1.35, **kw)


def band(ax, x0, x1, color, alpha=0.10, zorder=0, ymin=0, ymax=1):
    """纵向色带（用于高亮某一分类区间）"""
    from matplotlib.patches import Rectangle
    r = Rectangle((x0, ymin), x1 - x0, ymax - ymin, transform=ax.get_xaxis_transform(),
                  facecolor=color, edgecolor="none", alpha=alpha, zorder=zorder)
    ax.add_patch(r)
    return r


def halo_bar(ax, y, x0, x1, color, height=0.62, alpha_peak=0.30,
             n_strips=44, zorder=1):
    """
    "色彩晕"式误差表达（任务书🟢明确鼓励，替代传统误差棒）。
    在 [x0, x1] 区间内用密度渐变条带表达不确定度：
    中心 alpha 最高，向两端衰减到 ~0。纯呈现，不改数值。
    """
    xs = np.linspace(x0, x1, n_strips)
    mid = 0.5 * (x0 + x1)
    half = max(0.5 * (x1 - x0), 1e-12)
    for x in xs:
        t = abs(x - mid) / half            # 0(中心) → 1(端点)
        a = alpha_peak * (1.0 - t) ** 1.6
        w = (x1 - x0) / n_strips * 1.35
        ax.add_patch(plt.Rectangle((x - w / 2, y - height / 2), w, height,
                                   facecolor=color, edgecolor="none",
                                   alpha=max(a, 0.0), zorder=zorder))
    return ax


def gauss_halo(ax, y, center, sd, color, span=3.6, height=0.70,
               alpha_peak=0.30, n=180, zorder=1, x_clip=None):
    """
    高斯"色彩晕"：以 center±sd 画出零分布密度脊（alpha 随密度衰减）。
    用于 Fig 3 森林图表达 2000 次置换零分布，替代硬边误差棒。
    只用已交付的 null_mean / null_sd，不做任何新计算。
    """
    if not np.isfinite(sd) or sd <= 0:
        return None
    xs = np.linspace(center - span * sd, center + span * sd, n)
    dens = np.exp(-0.5 * ((xs - center) / sd) ** 2)
    dens = dens / dens.max()
    if x_clip is not None:
        lo, hi = x_clip
        keep = (xs >= lo) & (xs <= hi)
        xs, dens = xs[keep], dens[keep]
        if len(xs) < 3:
            return None
    for i in range(len(xs) - 1):
        a = alpha_peak * 0.5 * (dens[i] + dens[i + 1])
        ax.add_patch(plt.Rectangle(
            (xs[i], y - height * dens[i] / 2),
            xs[i + 1] - xs[i], height * dens[i],
            facecolor=color, edgecolor="none", alpha=float(a), zorder=zorder))
    return ax


def pfmt(p) -> str:
    """p 值格式化（纯呈现，不改数值）。经验 p 最小值 = 1/(1+nperm)。"""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "n/a"
    p = float(p)
    if p >= 0.001:
        return f"p = {p:.3f}".replace("0.", "0.")
    return f"p = {p:.1e}"


def er_fmt(er) -> str:
    if er is None:
        return "n/a"
    try:
        e = float(er)
    except (TypeError, ValueError):
        return str(er)
    if not np.isfinite(e):
        return "not computable"
    return f"{e:.2f}"


def sig_mark(er, p, er_min=2.0, p_max=0.01) -> bool:
    """round8/round9 criteria C3 预注册判定：ER>2 且 p<0.01 → 显著富集。"""
    try:
        return (float(er) > er_min) and (float(p) < p_max)
    except (TypeError, ValueError):
        return False


# ======================================================================
# 5. 三格式保存（PNG @dpi + SVG + PDF）
# ======================================================================
def save_figure(fig, outdir: Path, name: str, kind: str = "main",
                png_only: bool = False, dpi_override: int | None = None) -> dict:
    """
    同时输出 PNG + SVG + PDF（绘图指令通用规格）。
    kind='main' → 600 dpi；kind='line' → 1200 dpi（线稿/文字图）。
    返回 {ext: path}。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dpi = dpi_override or (DPI_MAIN if kind == "main" else DPI_LINE)

    paths = {}
    p_png = outdir / f"{name}.png"
    fig.savefig(p_png, dpi=dpi, facecolor="white")
    paths["png"] = p_png
    if not png_only:
        for ext in ("svg", "pdf"):
            p = outdir / f"{name}.{ext}"
            fig.savefig(p, facecolor="white")   # 矢量：dpi 无关
            paths[ext] = p
    plt.close(fig)
    return paths


# ======================================================================
# 6. source_data.tsv 读写器（脚本内不写死数字的唯一保障）
# ======================================================================
SD_COLS = ["block", "entity", "sub", "field", "value",
           "num", "den", "text", "provenance"]


class SourceTable:
    """
    每张图一份 source_data.tsv：长表（tidy）格式，列固定为 SD_COLS。
      block      —— 面板/数据块标识（A_glyph, B_rate, annot, ...）
      entity     —— 主体（目名、菌株 accession、家族名、模型名…）
      sub        —— 次级维度（数据集、层、指标…）
      field      —— 字段名
      value      —— 数值（字符串形式，绘图端 pd.to_numeric 解析）
      num / den  —— 比率类数值的分子/分母（保证"图上每个比例都能还原成计数"）
      text       —— 文字标注
      provenance —— 该行的上游文件 + 定位方式（可溯源，红线 R2）
    """

    def __init__(self, figure: str):
        self.figure = figure
        self.rows: list[dict] = []

    def add(self, block, entity="", sub="", field="", value="",
            num="", den="", text="", prov="") -> "SourceTable":
        self.rows.append(dict(block=block, entity=entity, sub=sub, field=field,
                              value=value, num=num, den=den, text=text,
                              provenance=prov))
        return self

    def add_num(self, block, entity, field, value, sub="", num="", den="",
                text="", prov="") -> "SourceTable":
        return self.add(block, entity, sub, field,
                        "" if value is None or (isinstance(value, float) and not np.isfinite(value)) else repr(float(value)) if isinstance(value, float) else str(value),
                        num, den, text, prov)

    def extend_df(self, block: str, df: pd.DataFrame, prov: str,
                  entity_col: str, field_map: dict | None = None,
                  sub_col: str | None = None, skip_cols=()) -> "SourceTable":
        """把一个 DataFrame 整体落进 source_data（逐格一行，零丢失）。"""
        fmap = field_map or {}
        for _, r in df.iterrows():
            ent = str(r[entity_col])
            sub = "" if sub_col is None else str(r[sub_col])
            for c in df.columns:
                if c in skip_cols or c == entity_col or (sub_col and c == sub_col):
                    continue
                v = r[c]
                if isinstance(v, float) and not np.isfinite(v):
                    sval = ""
                elif pd.isna(v):
                    sval = ""
                else:
                    sval = str(v)
                self.add(block, ent, sub, fmap.get(c, c), sval, prov=prov)
        return self

    def write(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(self.rows, columns=SD_COLS)
        header = (f"# source_data for {self.figure} | figure-round (presentation only)\n"
                  f"# rows={len(df)} | all plotted values are read back from this file by the plot stage\n"
                  f"# provenance column gives the upstream read-only file + locator for every value\n")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(header)
            df.to_csv(fh, sep="\t", index=False)
        return path

    def __len__(self):
        return len(self.rows)


class SourceFrame:
    """PLOT 阶段的只读视图：只认 source_data.tsv。"""

    def __init__(self, df: pd.DataFrame, path: Path):
        self.df = df
        self.path = path

    def block(self, name: str) -> pd.DataFrame:
        return self.df[self.df["block"] == name].reset_index(drop=True)

    def blocks(self) -> list[str]:
        return list(dict.fromkeys(self.df["block"].tolist()))

    @staticmethod
    def _num(s):
        return pd.to_numeric(s.replace({"": None, "nan": None}), errors="coerce")

    def values(self, block: str, field: str = "", entity: str | None = None,
               sub: str | None = None) -> pd.Series:
        d = self.block(block)
        if field:
            d = d[d["field"] == field]
        if entity is not None:
            d = d[d["entity"] == entity]
        if sub is not None:
            d = d[d["sub"] == sub]
        return self._num(d["value"])

    def one(self, block: str, field: str, entity: str | None = None,
            sub: str | None = None, default=np.nan) -> float:
        v = self.values(block, field, entity, sub)
        return float(v.iloc[0]) if len(v) else default

    def one_text(self, block: str, field: str, entity: str | None = None,
                 sub: str | None = None, default: str = "") -> str:
        d = self.block(block)
        if field:
            d = d[d["field"] == field]
        if entity is not None:
            d = d[d["entity"] == entity]
        if sub is not None:
            d = d[d["sub"] == sub]
        return str(d["text"].iloc[0]) if len(d) else default

    def matrix(self, block: str, index: str = "entity",
               columns: str = "field", value: str = "value") -> pd.DataFrame:
        d = self.block(block).copy()
        d[value] = self._num(d[value]) if value == "value" else d[value]
        return d.pivot_table(index=index, columns=columns, values=value, aggfunc="first")


def read_source(path: Path) -> SourceFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    for c in SD_COLS:
        if c not in df.columns:
            df[c] = ""
    return SourceFrame(df[SD_COLS], Path(path))


# ======================================================================
# 7. 上游表读取小工具
# ======================================================================
def read_upstream(key: str, **kw) -> pd.DataFrame:
    """读上游只读 TSV（跳过 # 注释头）。key 见 SRC 字典。"""
    p = SRC[key]
    if not p.exists():
        raise FileNotFoundError(f"[shared_style] missing upstream input: {key} -> {p}")
    return pd.read_csv(p, sep="\t", comment="#", **kw)


def rel(p) -> str:
    """把绝对路径转成相对上游输入目录的短名（用于 provenance 列，避免绝对路径）。"""
    p = Path(p)
    try:
        return str(p.relative_to(INPUT))
    except ValueError:
        return p.name


# ======================================================================
# 8. 画布 / 布局
# ======================================================================
def new_figure(w_mm: float, h_mm: float):
    fig = plt.figure(figsize=figsize(w_mm, h_mm))
    return fig


def grid_axes(fig, nrows, ncols, left, bottom, right, top,
              wspace=0.30, hspace=0.42, wratios=None, hratios=None):
    """用 GridSpec 精确控制版面（留白由调用方给比例，制造呼吸感）。"""
    gs = fig.add_gridspec(nrows, ncols, left=left, bottom=bottom, right=right, top=top,
                          wspace=wspace, hspace=hspace,
                          width_ratios=wratios, height_ratios=hratios)
    return gs


def rounded_panel(fig, x0, y0, x1, y1, fc=None, ec=None, r=0.012, alpha=1.0,
                  zorder=0, lw=0.6):
    """圆角面板底（用于给某一口径/某一结论加"视觉容器"，强化分区）。"""
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       transform=fig.transFigure, clip_on=False,
                       facecolor=fc or PALETTE["panel"],
                       edgecolor=ec or PALETTE["panel_edge"],
                       linewidth=lw, alpha=alpha, zorder=zorder)
    fig.patches.append(p)
    return p


__all__ = [n for n in dir() if not n.startswith("_")]
