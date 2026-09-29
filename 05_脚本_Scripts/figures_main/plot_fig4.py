#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig4.py — Fig 4 家族多功能性 + 功效分析
=============================================
【任务书要求】
(A) 条形图，横轴=家族，纵轴=Activities in Family 数，颜色区分 核心/次级/侧链；
    标注 "GH3:36 / GH1:35 / GH2:18"
(B) 功效曲线，横轴=阴性样本数，纵轴=期望 κ，画 κ≥0.4 阈值线，
    标注 "bottleneck is predictor quality, not sample size"

【数据来源说明（重要）】
任务书指定 (A) 的数据为「学生 TableA_cazy_activities.xlsx」，该文件**未随附**。
同一项目的已定稿数值在随附的《绘图轮全景任务表.MD》§2 表中（CAZy V5-2 家族
Activities 登记数与 Characterized 数），本脚本从该表解析，逐行记录 provenance。
GH54 / GH79 在 §2 中标注"定位待核"、无 Activities 数值 → 按红线 R1 留空并以
斜纹"no curated data"条呈现，绝不插补。GH5/GH66 属面板外扩展候选，不进入本面板。

(B) 数据：outputs_round8/figures/fig3_source_data.tsv（M1_binary 解析期望 κ 曲线）、
outputs_round8/tables/t6_scenario_grid.tsv（κ≥0.4 的最小 n、κ 峰值位置）。
当前预测器 κ = −0.056 引自《绘图轮全景任务表.MD》§1.2/§5。

运行：python plot_fig4.py / --plot-only
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

FIG_NAME = "fig4_family_multifunctionality_and_power"
OUT_DIR = S.MAIN_DIR
SOURCE_DATA = OUT_DIR / f"{FIG_NAME}_source_data.tsv"
WIDTH_MM, HEIGHT_MM = S.MM["double"], 112.0
KIND = "main"

ROLE_COLOR = {"core": S.SEM["role_core"], "secondary": S.SEM["role_secondary"],
              "sidechain": S.SEM["role_sidechain"], "nodata": S.SEM["role_nodata"]}
PANEL_ORDER = ["GH1", "GH2", "GH3", "GH27", "GH36", "GH51", "GH54", "GH78", "GH79", "GH106"]


# ======================================================================
# BUILD
# ======================================================================
def _parse_handoff_activities():
    """从《绘图轮全景任务表.MD》§2 表解析 Activities in Family / Characterized / 角色。"""
    txt = Path(S.SRC["handoff_md"]).read_text(encoding="utf-8")
    rows = {}
    for line in txt.splitlines():
        line = line.rstrip("\r")
        if not re.match(r"^GH\d+\t", line):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        fam, ec, act, ch, role_txt = f[0], f[1], f[2], f[3], f[4]
        role = ("core" if "核心" in role_txt else
                "secondary" if "次级" in role_txt else
                "sidechain" if "侧链" in role_txt else "outside")
        act_v = float(act) if re.fullmatch(r"\d+", act.strip()) else None
        ch_v = float(ch) if re.fullmatch(r"\d+", ch.strip()) else None
        rows[fam] = dict(ec=ec, activities=act_v, characterized=ch_v, role=role,
                         role_text=role_txt)
    return rows


def build() -> S.SourceTable:
    st = S.SourceTable("Fig 4 — CAZy family multifunctionality and statistical power")
    PROV_MD = S.rel(S.SRC["handoff_md"])

    act = _parse_handoff_activities()
    for fam in PANEL_ORDER:
        d = act.get(fam)
        if d is None:
            continue
        if d["activities"] is None:
            # §2 未给数值（GH54/GH79 定位待核）→ 留空 + 原因（红线 R1）
            st.add("A_activities", fam, d["role"] if d["role"] != "outside" else "nodata",
                   "activities", "", text="no curated activity count in handoff §2 (positioning pending)",
                   prov=PROV_MD)
        else:
            st.add("A_activities", fam, d["role"], "activities", str(d["activities"]), prov=PROV_MD)
        if d["characterized"] is not None:
            st.add("A_activities", fam, d["role"], "characterized", str(d["characterized"]), prov=PROV_MD)
        st.add("A_activities", fam, d["role"], "has_ec32121",
               "1" if "H15" in d["ec"] or "✅" in d["ec"] else "0", text=d["ec"], prov=PROV_MD)

    # --- B_power：期望 κ 曲线 -----------------------------------------
    pc = S.read_upstream("power_curve")
    pc = pc[pc["model"] == "M1_binary"]
    for _, r in pc.iterrows():
        st.add("B_power", f"Se{r['Se']}_Sp{r['Sp']}", "", "kappa", str(r["kappa"]),
               prov=S.rel(S.SRC["power_curve"]))
        st.add("B_power_n", f"Se{r['Se']}_Sp{r['Sp']}", str(r["n_neg"]), "kappa", str(r["kappa"]),
               prov=S.rel(S.SRC["power_curve"]))
    pg = S.read_upstream("power_grid")
    for _, r in pg.iterrows():
        if r["Se"] != r["Sp"]:
            continue
        st.add("B_grid", f"Se{r['Se']}", "", "min_n_kappa_ge_04", str(r["min_n_for_kappa_ge_04"]),
               prov=S.rel(S.SRC["power_grid"]))
        st.add("B_grid", f"Se{r['Se']}", "", "kappa_max", str(r["kappa_max"]), prov=S.rel(S.SRC["power_grid"]))
        st.add("B_grid", f"Se{r['Se']}", "", "n_at_kappa_max", str(r["n_at_kappa_max"]), prov=S.rel(S.SRC["power_grid"]))

    # 当前预测器 κ（引自交接档 §1.2/§5，非本图计算）
    txt = Path(S.SRC["handoff_md"]).read_text(encoding="utf-8")
    m = re.search(r"kappa\s*=\s*[−-]\s*0\.056", txt) or re.search(r"[−-]0\.056", txt)
    cur = -0.056 if m else None
    st.add("meta", "current_predictor_kappa", "", "value", str(cur),
           text="four-class strain-level predictor, chance level (handoff §1.2/§5)", prov=PROV_MD)
    st.add("meta", "kappa_target", "", "value", "0.4", prov="task spec (kappa >= 0.4 threshold line)")
    for i, f in enumerate(PANEL_ORDER):
        st.add("meta", "family_order", str(i), "name", text=f, prov="10-family panel (handoff §2)")

    st.write(SOURCE_DATA)
    print(f"[fig4] BUILD  -> {SOURCE_DATA.name}  ({len(st)} rows)")
    return st


# ======================================================================
# PLOT
# ======================================================================
def plot():
    S.apply_rc()
    sd = S.read_source(SOURCE_DATA)
    A = sd.block("A_activities").copy()
    fams = [t for t in (sd.one_text("meta", "name", entity="family_order", sub=str(i))
                        for i in range(len(PANEL_ORDER))) if t]
    kappa_target = sd.one("meta", "value", entity="kappa_target")
    cur_kappa = sd.one("meta", "value", entity="current_predictor_kappa")

    fig = S.new_figure(WIDTH_MM, HEIGHT_MM)
    gsA = fig.add_gridspec(1, 1, left=0.075, right=0.430, top=0.880, bottom=0.235)
    axA = fig.add_subplot(gsA[0, 0])
    gsB = fig.add_gridspec(1, 1, left=0.530, right=0.985, top=0.880, bottom=0.235)
    axB = fig.add_subplot(gsB[0, 0])

    # ---------------- Panel A：Activities in Family ----------------
    y = np.arange(len(fams))[::-1]
    for i, fam in enumerate(fams):
        r = A[(A["entity"] == fam) & (A["field"] == "activities")]
        role = r["sub"].iloc[0] if len(r) else "nodata"
        col = ROLE_COLOR.get(role, S.SEM["role_nodata"])
        val = pd.to_numeric(r["value"], errors="coerce")
        ch = A[(A["entity"] == fam) & (A["field"] == "characterized")]
        chv = pd.to_numeric(ch["value"], errors="coerce")
        chv = float(chv.iloc[0]) if len(chv) and np.isfinite(chv.iloc[0]) else None

        if len(val) and np.isfinite(val.iloc[0]):
            v = float(val.iloc[0])
            axA.barh(y[i], v, height=0.62, color=col, alpha=0.90,
                     edgecolor="white", linewidth=0.5, zorder=3)
            lbl = f"{v:g}"
            bold = fam in ("GH1", "GH2", "GH3")
            axA.text(v + 0.6, y[i], lbl, va="center", ha="left",
                     fontsize=S.FS["note"], color=S.PALETTE["ink"],
                     fontweight="bold" if bold else "normal")
            if chv is not None:
                axA.text(48.6, y[i], f"char. n = {chv:g}", va="center", ha="right",
                         fontsize=S.FS["note"], color=S.PALETTE["ink_faint"])
        else:
            # 数据缺失：斜纹留空条 + 原因（不插补）
            axA.barh(y[i], 38.0, height=0.62, color="none", edgecolor=S.PALETTE["ink_faint"],
                     linewidth=0.6, hatch="////", zorder=3)
            axA.text(0.6, y[i], "no curated activity count", va="center", ha="left",
                     fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], style="italic")
        axA.text(-0.6, y[i], fam, transform=axA.get_yaxis_transform(), ha="right", va="center",
                 fontsize=S.FS["tick"], color=col if role != "nodata" else S.PALETTE["ink_faint"],
                 fontweight="bold" if fam in ("GH1", "GH2", "GH3") else "normal")

    axA.axvline(38.5, color=S.PALETTE["grid"], lw=0.6, zorder=1)
    axA.set_ylim(-0.7, len(fams) - 0.3)
    axA.set_xlim(0, 49)
    axA.set_xticks([0, 10, 20, 30, 40])
    axA.set_yticks([])
    axA.set_xlabel("distinct activities registered in CAZy V5-2 per family", fontsize=S.FS["axis"])
    S.style_ax(axA, grid="x", keep_spines=("bottom",))
    S.panel_label(axA, "A", "Family multifunctionality   (GH3: 36 · GH1: 35 · GH2: 18)")
    # 角色图例
    from matplotlib.patches import Patch
    axA.legend(handles=[Patch(fc=ROLE_COLOR["core"], label="core (EC 3.2.1.21)"),
                        Patch(fc=ROLE_COLOR["secondary"], label="secondary"),
                        Patch(fc=ROLE_COLOR["sidechain"], label="side-chain"),
                        Patch(fc="none", ec=S.PALETTE["ink_faint"], hatch="////", label="no curated data")],
               loc="center left", fontsize=S.FS["note"], ncol=1, columnspacing=0.9,
               handletextpad=0.4, labelspacing=0.30, bbox_to_anchor=(0.33, 0.60))

    # ---------------- Panel B：power curves ----------------
    B = sd.block("B_power_n")
    curves = sorted({e for e in B["entity"]})
    # 只画 Se==Sp 的对角线（与 round8 fig3 同口径），按 Se/Sp 用单色相明度梯度
    diag = [c for c in curves if c.replace("Se", "").split("_Sp")[0] == c.split("_Sp")[1]]
    ses = sorted({float(c.split("_Sp")[1]) for c in diag}, reverse=True)
    ramp = plt.get_cmap("Blues_r")
    for k, se in enumerate(ses):
        key = f"Se{se:g}_Sp{se:g}"
        d = B[B["entity"] == key].copy()
        nn = pd.to_numeric(d["sub"], errors="coerce")
        d = d.assign(nn=nn).sort_values("nn")
        kap = pd.to_numeric(d["value"], errors="coerce")
        is_chance = abs(se - 0.5) < 1e-9
        col = S.PALETTE["grey"] if is_chance else ramp(0.15 + 0.75 * (k / max(len(ses) - 2, 1)))
        axB.plot(d["nn"].clip(lower=1), kap,
                 linestyle=":" if is_chance else "-",
                 color=col, lw=1.5 if not is_chance else 1.0,
                 label=None, zorder=3)
        axB.text(d["nn"].iloc[-1], kap.iloc[-1], f"{se:g}", fontsize=S.FS["note"],
                 color=col, va="center", ha="left",
                 fontweight="bold" if not is_chance else "normal")

    axB.axhline(kappa_target, color=S.PALETTE["ink_soft"], lw=0.9, ls="--", zorder=2)
    axB.text(1.15, kappa_target + 0.02, "target κ = 0.4", fontsize=S.FS["note"],
             color=S.PALETTE["ink_soft"], va="bottom")
    if np.isfinite(cur_kappa):
        axB.axhline(cur_kappa, color=S.SEM["not_detected"], lw=0.9, ls=":", zorder=2)
        axB.text(1.15, cur_kappa - 0.03, f"current family-level predictor κ = {cur_kappa:.3f} (chance)",
                 fontsize=S.FS["note"], color=S.SEM["not_detected"], va="top")

    # κ≥0.4 的最小 n（已交付）标在阈值线上
    G = sd.block("B_grid")
    _minn_parts = []
    gmin = G[G["field"] == "min_n_kappa_ge_04"].set_index("entity")["value"]
    for ent, v in gmin.items():
        se = float(ent.replace("Se", ""))
        if v == "never" or not str(v).replace(".", "", 1).isdigit():
            continue
        n = float(v)
        col = ramp(0.15 + 0.75 * ((ses.index(se) if se in ses else 0) / max(len(ses) - 2, 1)))
        axB.scatter([n], [kappa_target], s=22, color=col, edgecolors=S.PALETTE["ink"],
                    linewidths=0.5, zorder=5, marker="v")
        _minn_parts.append(f"{se:g}\u2192{n:g}")

    # κ 峰值（非单调）注记
    gmax = G[G["field"] == "n_at_kappa_max"].set_index("entity")["value"]
    if "Se0.9" in gmax.index:
        nmax = float(gmax.loc["Se0.9"])
        axB.axvline(nmax, color=S.PALETTE["ink_faint"], lw=0.7, ls=(0, (2, 2)), zorder=1)


    _never = [e.replace("Se", "") for e, v in gmin.items() if str(v) == "never"]
    axB.set_xscale("log")
    axB.set_xlim(1, 22000)
    axB.set_ylim(-0.14, 1.02)
    axB.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axB.set_xlabel("confirmed-negative samples added (positives fixed at 952)", fontsize=S.FS["axis"])
    axB.set_ylabel("expected Cohen κ", fontsize=S.FS["axis"])
    S.style_ax(axB, grid="y")
    S.panel_label(axB, "B", "Power: κ vs confirmed-negative sample size (Se = Sp curves labelled)")

    # ---------------- 强制注记 ----------------
    gsN = fig.add_gridspec(1, 1, left=0.075, right=0.985, top=0.150, bottom=0.012)
    axN = fig.add_subplot(gsN[0, 0]); axN.axis("off")
    axN.text(0, 1.0,
        "Bottleneck is predictor quality, not sample size: at the current 952:2 ascertainment bias a predictor at chance level never reaches κ = 0.4 at any n, "
        "whereas Se = Sp ≥ 0.75 needs only ~10² confirmed negatives.  κ(n) is non-monotonic (it peaks near balanced prevalence).  "
        "Delivered minimum n for E[κ] ≥ 0.4 (Se=Sp): " + ", ".join(_minn_parts)
        + ("; " + "/".join(_never) + " never reach 0.4." if _never else "."),
        fontsize=S.FS["note"], color=S.PALETTE["ink"], va="top", fontweight="bold")
    axN.text(0, 0.52,
        "Panel A: a family that registers 4–36 distinct activities cannot resolve substrate specificity from copy number alone — this is the structural reason "
        "family-level genomic prediction of Rb1→CK conversion fails (main line B).  Activities/Characterized counts: CAZy V5-2 as tabulated in the project handoff §2 "
        "(the student TableA_cazy_activities.xlsx was not delivered; GH54/GH79 have no curated count and are left blank, not imputed).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_soft"], va="top")
    axN.text(0, 0.10,
        "Panel B curves: analytic expectation of κ under a binary M1 model (outputs_round8/figures/fig3_source_data.tsv); ▼ = delivered minimum n for E[κ] ≥ 0.4 (t6_scenario_grid.tsv).",
        fontsize=S.FS["note"], color=S.PALETTE["ink_faint"], va="top")

    S.figure_title(fig, "Why family-level copy number cannot predict function: multifunctionality and the power limit",
                   x=0.075, y=0.988)

    paths = S.save_figure(fig, OUT_DIR, FIG_NAME, kind=KIND)
    print("[fig4] PLOT   ->", ", ".join(sorted(paths)))
    return paths


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        build()
    plot()
