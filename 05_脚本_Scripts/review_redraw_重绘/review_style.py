#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
review_style.py — 审稿补数据图（Fig R1/R2/R3）统一绘图样式
=================================================================
设计原则（全三图一致）：
  · 无衬线字体：优先 Arial，回退 Liberation Sans / DejaVu Sans；
  · 配色 Okabe–Ito 色盲安全色板（与原论文绘图轮 shared_style.py 同源）；
  · 同一模块跨图同色：GH3=绿、GH5=紫粉；model P=蓝、model C=橙；
    表达证据四级分类 = 蓝/天蓝/橙/灰；
  · 无图内图注：图面只保留轴标签、刻度、短图例与面板标签 A/B/C；
  · 输出 SVG + PDF + PNG(300 dpi)，双栏幅宽 183 mm。
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------- 字体
for _f in ("LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf",
           "LiberationSans-Italic.ttf", "LiberationSans-BoldItalic.ttf"):
    _p = Path.home() / ".fonts" / _f
    if _p.exists():
        mpl.font_manager.fontManager.addfont(str(_p))

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 7.0,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7.0,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.6,
    "grid.linewidth": 0.4,
    "axes.grid": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",          # SVG 保留可编辑文本
    "pdf.fonttype": 42,              # PDF 嵌入 TrueType
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.01,
})

# ---------------------------------------------------------------- 色板（Okabe–Ito，跨图一致）
C = {
    "ink":       "#1A1A1A",
    "ink_soft":  "#5C5C5C",
    "ink_faint": "#8C8C8C",
    "grid":      "#E8E8E8",
    "panel":     "#F6F6F6",
    # 家族（全套图一致）
    "GH3":       "#009E73",   # green
    "GH5":       "#CC79A7",   # reddish purple
    "GH3_l":     "#B4E0D2",
    "GH5_l":     "#EED3E2",
    # 零模型（全套图一致）
    "modelP":    "#0072B2",   # blue
    "modelC":    "#E69F00",   # orange
    "modelP_l":  "#B8D9EA",
    "modelC_l":  "#F7DFAE",
    # 宏转录组表达四级分类（全套图一致）
    "expr_exact":   "#0072B2",   # near-exact ≥95 %
    "expr_high":    "#56B4E9",   # high-identity ≥90 %
    "expr_homolog": "#E69F00",   # homolog 60–90 %
    "expr_nohit":   "#B8B8B8",   # no hit
    # 强调
    "vermillion": "#D55E00",
    "blue":       "#0072B2",
    "grey":       "#646464",
}

MM = 1 / 25.4  # mm → inch


def new_fig(width_mm: float, height_mm: float):
    return plt.figure(figsize=(width_mm * MM, height_mm * MM))


def panel_label(ax, label: str, x: float = -0.16, y: float = 1.06):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=8.5, fontweight="bold", va="top", ha="left", color=C["ink"])


def strip_axes(ax, left=True, bottom=True):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


def save_fig(fig, outdir: Path, name: str) -> list[Path]:
    """统一三格式输出：SVG + PDF + PNG(300 dpi)。"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("svg", "pdf", "png"):
        p = outdir / f"{name}.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None)
        paths.append(p)
    plt.close(fig)
    return paths


def short_wp(wp: str) -> str:
    """WP_010801083.1 → …010801083"""
    core = wp.replace("WP_", "").split(".")[0]
    return "…" + core
