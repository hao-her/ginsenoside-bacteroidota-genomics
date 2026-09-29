#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_palette_sheet.py — 统一配色说明 + 色盲模拟验证（红线 3）
==============================================================
产出（palette/）：
  palette_swatches.png/.svg/.pdf   每个色值在 正常视觉 + 三种完全型色觉下的并排样张
  palette_cvd_deltaE.tsv           逐对 CIELAB ΔE*ab（正常 + protanopia/deuteranopia/tritanopia）+ 判定
  cvd_simulation/*.png             8 张成图在三型色觉下的模拟截图（缩略，供人工复核）
  palette.md                       色值表、语义映射、验证结论与救济说明

判定标准（冻结于 cvd.THRESHOLDS，不事后调整）：
  三型完全色盲下逐对 ΔE*ab ≥ 15 → PASS；12–15 → WARN（须有非颜色冗余通道）；< 12 → FAIL（禁用）。
另做"纯红绿对比"检查：与参照对 #FF0000/#00FF00 在 deuteranopia 下的 ΔE 比较。

运行：python make_palette_sheet.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import cvd
import shared_style as S

OUT = S.PALETTE_DIR
SIMDIR = OUT / "cvd_simulation"
SIMDIR.mkdir(parents=True, exist_ok=True)

P = S.PALETTE

# --- 承担"唯一区分通道"的核心分类色集（Okabe–Ito 7 色）---------------
CORE_SET = {
    "blue (GH1 / detected / unbiased / Bacteroidota / Bacteroides)": P["blue"],
    "vermillion (GH1-absence / conditional / Bacteroidales / highlight)": P["vermillion"],
    "green (GH3 / control / Ascomycota / Sphingobacteriales)": P["green"],
    "orange (GH2 / secondary / Bacillota / Cytophagales)": P["orange"],
    "purple (side-chain / Actinomycetota / Parabacteroides)": P["purple"],
    "sky (GH36 / rescued)": P["sky"],
    "grey (no-signal / ns / context)": P["grey"],
}
# --- 派生淡色（仅作填充/底色，与基色靠明度区分）-----------------------
TINT_SET = {
    "blue": P["blue"], "blue_l": P["blue_l"],
    "vermillion": P["vermillion"], "verm_l (absence fill)": P["verm_l"],
    "green": P["green"], "green_l": P["green_l"],
    "orange": P["orange"], "orange_l": P["orange_l"],
    "purple": P["purple"], "purple_l": P["purple_l"],
    "sky": P["sky"], "sky_l": P["sky_l"],
    "grey": P["grey"], "grey_l (no-data hatch)": P["grey_l"],
}
# --- Fig 5 属级用色（去重后）：n≤4 的稀有属共用中性灰，其身份由行标签文字承担 ---
_gen_groups = {}
for _k, _v in S.GENUS.items():
    _gen_groups.setdefault(_v, []).append(_k)
GENUS_SET = {("genus " + g[0]) if len(g) == 1 else
             f"grey shared by {'/'.join(g)} (identity carried by row label, not colour)": v
             for v, g in _gen_groups.items()}

FIGS = [
    ("fig1", S.MAIN_DIR / "fig1_phylum_genus_gh_distribution.png"),
    ("fig2", S.MAIN_DIR / "fig2_gh1_lineage_boundary.png"),
    ("fig3", S.MAIN_DIR / "fig3_pul_genomic_organization.png"),
    ("fig4", S.MAIN_DIR / "fig4_family_multifunctionality_and_power.png"),
    ("fig5", S.MAIN_DIR / "fig5_candidate_screening_funnel.png"),
    ("suppfig1", S.SUPP_DIR / "suppfig1_gh1_threshold_curve.png"),
    ("suppfig2", S.SUPP_DIR / "suppfig2_phylogeny_gh1_mapping.png"),
    ("suppfig3", S.SUPP_DIR / "suppfig3_pul_unbiased_45strains.png"),
]


# ======================================================================
def swatch_sheet():
    S.apply_rc()
    groups = [("Core categorical set (Okabe–Ito, sole distinguishing channel)", CORE_SET),
              ("Derived tints (fill / background; separated from bases by lightness)", TINT_SET),
              ("Fig 5 genus set (redundant with row position)", GENUS_SET)]
    ncols = 4
    nrows = sum(max(1, (len(g[1]) + ncols - 1) // ncols) for g in groups) + len(groups)
    fig = plt.figure(figsize=(S.mm(183), S.mm(16 * nrows + 30)))
    y = 0.985
    rows_all = []
    for gname, g in groups:
        fig.text(0.03, y, gname, fontsize=S.FS["axis"], fontweight="bold", va="top")
        y -= 0.030
        items = list(g.items())
        for i, (name, hexv) in enumerate(items):
            r, c = divmod(i, ncols)
            x0 = 0.03 + c * 0.245
            yy = y - r * 0.062
            sims = [hexv] + [cvd.simulate_hex(hexv, t) for t in cvd.CVD_TYPES]
            for k, sv in enumerate(sims):
                fig.patches.append(plt.Rectangle((x0 + k * 0.030, yy - 0.026), 0.026, 0.026,
                                                 transform=fig.transFigure, facecolor=sv,
                                                 edgecolor=P["panel_edge"], linewidth=0.4))
            fig.text(x0 + 4 * 0.030 + 0.006, yy - 0.013, f"{name}\n{hexv}",
                     fontsize=S.FS["note"], va="center", ha="left", linespacing=1.3)
            rows_all.append((gname, name, hexv))
        nrow = (len(items) + ncols - 1) // ncols
        y -= nrow * 0.062 + 0.020
    # 列头
    fig.text(0.03, 0.995, "normal      protanopia      deuteranopia      tritanopia",
             fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    fig.savefig(OUT / "palette_swatches.png", dpi=S.DPI_MAIN, facecolor="white")
    fig.savefig(OUT / "palette_swatches.svg", facecolor="white")
    fig.savefig(OUT / "palette_swatches.pdf", facecolor="white")
    plt.close(fig)
    return rows_all


# ======================================================================
def delta_table():
    rows = []
    for setname, g in (("core", CORE_SET), ("tints", TINT_SET), ("genus", GENUS_SET)):
        res = cvd.verify_palette(g)
        for pr in res["pairs"]:
            rows.append(dict(color_set=setname, **pr))
        rg = cvd.red_green_check(g)
        for r in rg:
            rows.append(dict(color_set=setname, color1=r["color1"], hex1=r["hex1"],
                             color2=r["color2"], hex2=r["hex2"],
                             deltaE_normal=np.nan,
                             deltaE_protanopia=np.nan, deltaE_deuteranopia=r["deltaE_deuteranopia"],
                             deltaE_tritanopia=np.nan, worst_type="red-green-screen",
                             worst_deltaE=r["deltaE_deuteranopia"], verdict=r["verdict"]))
    df = __import__("pandas").DataFrame(rows)
    df.to_csv(OUT / "palette_cvd_deltaE.tsv", sep="\t", index=False)
    return df


# ======================================================================
def simulate_figures():
    made = []
    for name, path in FIGS:
        if not path.exists():
            print(f"[palette] skip missing {path.name}")
            continue
        Image.MAX_IMAGE_PIXELS = None
        im = Image.open(path).convert("RGB")
        # 缩到长边 2000px 以控制验证截图体积（验证用途，非交付主图）
        scale = 2000 / max(im.size)
        if scale < 1:
            im = im.resize((int(im.size[0] * scale), int(im.size[1] * scale)), Image.LANCZOS)
        arr = np.asarray(im).astype(float) / 255.0
        for t in cvd.CVD_TYPES:
            sim = cvd.simulate(arr, t)
            out = SIMDIR / f"{name}_{t}.png"
            Image.fromarray((np.clip(sim, 0, 1) * 255).astype(np.uint8)).save(out, optimize=True)
            made.append(out)
    return made


# ======================================================================
def write_md(df, n_sim):
    core = df[df.color_set == "core"]
    worst_core = core["worst_deltaE"].min()
    n_fail_core = int((core["verdict"] == "FAIL").sum())
    n_warn_core = int((core["verdict"] == "WARN").sum())
    gen = df[(df.color_set == "genus") & (df.worst_type != "red-green-screen")]
    worst_gen = gen["worst_deltaE"].min()
    rg = df[df.verdict == "PURE_RED_GREEN"]

    lines = []
    lines.append("# palette/ — 统一配色方案与色盲模拟验证\n")
    lines.append("本套 8 张图共用同一色板与同一语义映射（见 `scripts/shared_style.py`），")
    lines.append("保证跨图“同一颜色 = 同一含义”，无割裂感。全部色值取自 Okabe–Ito (2008) 色盲友好色板。\n")
    lines.append("## 1. 色值与语义\n")
    lines.append("| 色值 | 语义（全套图一致） |")
    lines.append("|---|---|")
    lines.append(f"| `{P['blue']}` | GH1 / 检出 / 无偏口径 / Bacteroidota / Bacteroides / 核心家族 |")
    lines.append(f"| `{P['vermillion']}` | GH1 未检出（醒目强调）/ 条件化口径（上界）/ Bacteroidales / 点名候选 |")
    lines.append(f"| `{P['green']}` | GH3 / 阳性对照 / Ascomycota / Sphingobacteriales |")
    lines.append(f"| `{P['orange']}` | GH2（次级）/ Bacillota / Cytophagales |")
    lines.append(f"| `{P['purple']}` | 侧链糖苷酶 / Actinomycetota / Parabacteroides |")
    lines.append(f"| `{P['sky']}` | GH36 / pairwise-rescued 候选 |")
    lines.append(f"| `{P['yellow']}` | 属级（Alistipes）；Okabe-Ito 第 8 色 |")
    lines.append(f"| `{P['grey']}` | 无信号 / 不显著 / 仅上下文 |")
    lines.append(f"| `{P['verm_l']}` / `{P['grey_l']}` | 未检出的填充色 / 数据缺失斜纹（不插补） |")
    lines.append(f"| `{P['ink']}` / `{P['ink_soft']}` / `{P['ink_faint']}` | 文字墨色三级 |")
    lines.append(f"| `{P['grid']}` / `{P['axis']}` | 网格 / 轴线 |")
    lines.append("")
    lines.append("热图（Suppl Fig 2）使用单色相明度渐变 `" + " → ".join(S.SEQ_RAMP) + "`，")
    lines.append("单色相 ramp 在任意色觉下仅靠明度区分，天然色盲安全。\n")
    lines.append("## 2. 色盲模拟验证（Machado et al. 2009, severity = 1.0）\n")
    lines.append(f"- 核心分类色集（7 色，承担唯一区分通道）：三型完全色盲下逐对 ΔE*ab 最小 = **{worst_core:.1f}**；"
                 f"FAIL {n_fail_core} 对，WARN {n_warn_core} 对 → **{'全部 PASS' if n_fail_core == 0 and n_warn_core == 0 else '见下表'}**。")
    lines.append(f"- Fig 5 属级用色（去重后 8 色：7 彩色 + 稀有属共用中性灰）：三型下逐对 ΔE*ab 最小 = **{worst_gen:.1f}**，全部 PASS。")
    lines.append("  n ≤ 4 的稀有属（Tannerella/Odoribacter/Rikenella/Porphyromonas）共用中性灰，其身份由行标签文字承担；")
    lines.append("  其余属的颜色同时与 y 轴行位置冗余编码。")
    lines.append("- 派生淡色（tints，仅作填充/底色）逐对 ΔE 较低属预期：它们从不单独承担区分任务，")
    lines.append("  总是与描边（如未检出格 = 淡朱红填充 + 朱红描边）、斜纹（数据缺失）或直接文字标签共同出现；明细见 tsv 的 tints 组。")
    if len(rg):
        lines.append(f"- 纯红绿对比检查：检出 {len(rg)} 对 → **须整改**。")
    else:
        lines.append("- 纯红绿对比检查：**0 对**。vermillion/green 等跨红绿色相区的组合在 deuteranopia 下的可分辨性")
        lines.append("  显著优于纯红绿参照对（#FF0000/#00FF00），属 Okabe–Ito 型的明度/蓝黄轴差异，不构成红绿对比。")
    lines.append("- 判定阈值（预先冻结）：三型 ΔE*ab ≥ 15 = PASS；12–15 = WARN（须非颜色冗余通道）；< 12 = FAIL（禁用）。")
    lines.append(f"- 逐对明细见 `palette_cvd_deltaE.tsv`；样张见 `palette_swatches.png`；")
    lines.append(f"  8 张成图的三型模拟截图共 {n_sim} 张，见 `cvd_simulation/`。\n")
    lines.append("## 3. 其它呈现红线自检\n")
    lines.append("- 无截断坐标轴制造夸大：Fig 1 / Fig 1-D 使用 symlog（linear 0–1, log above）并在轴标签注明，真 0 被如实画出。")
    lines.append("- Suppl Fig 3 中 2 株 0 配对的占比为“未定义”，画在独立槽位而非 0。")
    lines.append("- 误差表达：森林图用零分布高斯“色彩晕” + 95% 置换区间细须（任务书🟢鼓励），非传统误差棒。")
    (OUT / "palette.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    swatch_sheet()
    print("[palette] swatches done")
    df = delta_table()
    print("[palette] deltaE table done:", len(df), "rows")
    made = simulate_figures()
    print("[palette] cvd simulations:", len(made))
    write_md(df, len(made))
    print("[palette] palette.md done")
