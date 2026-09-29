#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cvd.py — 色觉障碍（CVD）模拟与配色可分辨性验证
================================================
红线 3：任何配色方案必须通过色盲模拟检验（deuteranopia / protanopia /
tritanopia 三型），禁止纯红绿对比。

实现
----
1) Machado, Oliveira & Fernandes (2009) 的线性变换矩阵，severity = 1.0
   （完全型 dichromacy，比 anomalous trichromacy 更严格 —— 通过完全型即通过轻度型）。
   变换在 **线性 RGB** 空间进行，故先做 sRGB gamma 解码，变换后再编码回 sRGB。
2) 可分辨性用 CIELAB ΔE*ab 量化：在 D65 白点下把原色与模拟色转 Lab，
   计算每一对颜色的 ΔE。判据（预先固定，见下 THRESHOLDS）。

判定标准（本文件冻结，不事后调整）
----------------------------------
- ΔE_pair_original  ≥ 20  ：原色下"清晰可辨"
- ΔE_pair_cvd       ≥ 15  ：三型完全色盲下仍"可辨"（PASS）
- 12 ≤ ΔE < 15            ：边际（WARN）→ 必须靠非颜色冗余通道救济
                            （直接文字标签 / 不同 marker / 位置分组），并在 README 注明
- ΔE < 12                 ：FAIL → 该两色不得在同一图中承担区分任务

注：数据色若同时用位置/形状/直接标签冗余编码，WARN 可接受（这是本套图的做法）。
"""

from __future__ import annotations

import numpy as np

# ----------------------------------------------------------------------
# Machado et al. 2009, severity = 1.0（线性 RGB 空间）
# ----------------------------------------------------------------------
CVD_MATRICES = {
    "protanopia": np.array([
        [0.152286, 1.052583, -0.204868],
        [0.114503, 0.786281,  0.099216],
        [-0.003882, -0.048116, 1.051998],
    ]),
    "deuteranopia": np.array([
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501,  0.047413],
        [-0.011820, 0.042940, 0.968881],
    ]),
    "tritanopia": np.array([
        [1.255528, -0.076749, -0.178779],
        [-0.078411, 0.930809,  0.147602],
        [0.004733, 0.691367,  0.303900],
    ]),
}
CVD_TYPES = ("protanopia", "deuteranopia", "tritanopia")

# sRGB (D65) → XYZ
RGB2XYZ = np.array([
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
])
# D65 白点
D65 = np.array([0.95047, 1.00000, 1.08883])

THRESHOLDS = {"pass": 15.0, "warn": 12.0, "clear_original": 20.0}


# ----------------------------------------------------------------------
# 颜色空间转换
# ----------------------------------------------------------------------
def hex2rgb(h: str) -> np.ndarray:
    h = h.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    return np.array([int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)])


def rgb2hex(rgb) -> str:
    a = np.clip(np.asarray(rgb, dtype=float), 0, 1)
    return "#{:02X}{:02X}{:02X}".format(*[int(round(v * 255)) for v in a])


def _srgb_to_linear(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(c: np.ndarray) -> np.ndarray:
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * np.clip(c, 0, None) ** (1 / 2.4) - 0.055)


def simulate(rgb, cvd_type: str) -> np.ndarray:
    """rgb: (..., 3) 0–1 sRGB。返回该色觉类型下"看起来"的 sRGB 颜色。"""
    if cvd_type == "normal":
        return np.asarray(rgb, dtype=float)
    M = CVD_MATRICES[cvd_type]
    arr = np.atleast_2d(np.asarray(rgb, dtype=float))
    lin = _srgb_to_linear(arr)
    out = lin @ M.T
    return np.clip(_linear_to_srgb(out), 0, 1)


def simulate_hex(hx: str, cvd_type: str) -> str:
    return rgb2hex(simulate(hex2rgb(hx), cvd_type)[0])


def rgb2lab(rgb) -> np.ndarray:
    """sRGB (0–1) → CIELAB (D65)。支持 (..., 3)。"""
    arr = np.atleast_2d(np.asarray(rgb, dtype=float))
    lin = _srgb_to_linear(np.clip(arr, 0, 1))
    xyz = lin @ RGB2XYZ.T
    xyz = xyz / D65

    def f(t):
        return np.where(t > 0.008856451679, np.cbrt(np.clip(t, 0, None)),
                        t / 0.12841854976 + 0.13793103448)

    fx, fy, fz = f(xyz[:, 0]), f(xyz[:, 1]), f(xyz[:, 2])
    L = 116 * fy - 16
    a = 500 * (fx - fy)
    b = 200 * (fy - fz)
    return np.column_stack([L, a, b])


def hex2lab(hx: str) -> np.ndarray:
    return rgb2lab(hex2rgb(hx))[0]


def delta_e(lab1, lab2) -> float:
    """CIE76 ΔE*ab（对"能否分辨"足够；非色差工程用途）。"""
    return float(np.sqrt(np.sum((np.asarray(lab1) - np.asarray(lab2)) ** 2)))


def luminance(rgb) -> float:
    """WCAG 相对亮度（0–1），用于文字/背景对比度检查。"""
    lin = _srgb_to_linear(np.clip(np.asarray(rgb, dtype=float), 0, 1))
    return float(np.dot([0.2126, 0.7152, 0.0722], lin))


def contrast_ratio(hex1: str, hex2: str) -> float:
    l1, l2 = luminance(hex2rgb(hex1)), luminance(hex2rgb(hex2))
    hi, lo = max(l1, l2), min(l1, l2)
    return (hi + 0.05) / (lo + 0.05)


# ----------------------------------------------------------------------
# 调色板级验证
# ----------------------------------------------------------------------
def verify_palette(colors: dict, labels: dict | None = None) -> dict:
    """
    colors: {name: '#RRGGBB'}
    返回 {'pairs': [ {...}, ... ], 'summary': {...}}
    每对给出原色 ΔE 与三型 CVD 下的 ΔE，并给 verdict。
    """
    names = list(colors)
    labs_normal = {n: hex2lab(colors[n]) for n in names}
    labs_cvd = {t: {n: rgb2lab(simulate(hex2rgb(colors[n]), t)[0])[0]
                    for n in names} for t in CVD_TYPES}
    sim_hex = {t: {n: simulate_hex(colors[n], t) for n in names} for t in CVD_TYPES}

    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            de0 = delta_e(labs_normal[n1], labs_normal[n2])
            dec = {t: delta_e(labs_cvd[t][n1], labs_cvd[t][n2]) for t in CVD_TYPES}
            worst_type = min(dec, key=dec.get)
            worst = dec[worst_type]
            if worst >= THRESHOLDS["pass"]:
                verdict = "PASS"
            elif worst >= THRESHOLDS["warn"]:
                verdict = "WARN"
            else:
                verdict = "FAIL"
            pairs.append(dict(
                color1=n1, hex1=colors[n1], color2=n2, hex2=colors[n2],
                deltaE_normal=round(de0, 2),
                **{f"deltaE_{t}": round(dec[t], 2) for t in CVD_TYPES},
                worst_type=worst_type, worst_deltaE=round(worst, 2),
                verdict=verdict,
            ))

    n_pass = sum(p["verdict"] == "PASS" for p in pairs)
    n_warn = sum(p["verdict"] == "WARN" for p in pairs)
    n_fail = sum(p["verdict"] == "FAIL" for p in pairs)
    return {
        "pairs": pairs,
        "simulated_hex": sim_hex,
        "summary": dict(n_pairs=len(pairs), n_pass=n_pass, n_warn=n_warn, n_fail=n_fail,
                        min_deltaE_cvd=min((p["worst_deltaE"] for p in pairs), default=float("nan")),
                        all_pass=(n_fail == 0)),
    }


def red_green_check(colors: dict) -> list[dict]:
    """
    红线检查：禁止"纯红绿对比"。

    仅靠色相分区判定会误伤 Okabe–Ito —— 该色板本就是为替代红绿对比而设计
    （#D55E00 / #009E73 是任务书明确允许的色值）。因此这里用**可分辨性**判定：

      1) 先按 CIELCh 色相角找出"红区 vs 绿区"的候选对；
      2) 再算该对在 deuteranopia 下的 ΔE*ab，与"纯红绿"参照对
         #FF0000 / #00FF00 在同型模拟下的 ΔE 比较（该参照值在下方 reference 键给出）。

    verdict:
      PURE_RED_GREEN      —— 德型 ΔE < 参照值：等同纯红绿对比，禁止（红线 3）
      OK_HUE_SHIFTED      —— 德型 ΔE ≥ 参照值：虽跨红/绿色相区，但在全色盲下
                             仍可分辨（Okabe–Ito 型明度/蓝黄轴差异），允许使用
    """
    ref_a = simulate(hex2rgb("#FF0000"), "deuteranopia")[0]
    ref_b = simulate(hex2rgb("#00FF00"), "deuteranopia")[0]
    ref_de = delta_e(rgb2lab(ref_a)[0], rgb2lab(ref_b)[0])

    out = []
    names = list(colors)
    lab = {n: hex2lab(colors[n]) for n in names}
    lab_de = {n: rgb2lab(simulate(hex2rgb(colors[n]), "deuteranopia")[0])[0] for n in names}

    def lch(L):
        C = np.hypot(L[1], L[2])
        h = np.degrees(np.arctan2(L[2], L[1])) % 360
        return C, h

    def zone(h):
        # 仅"接近纯红/纯绿"的色相才进入筛查；Okabe-Ito 的 vermillion(≈56°) 与
        # bluish-green(≈164°)、reddish-purple(≈344°) 均不在此范围 —— 它们正是为
        # 替代红绿对比而设计的色值。
        if (h <= 20) or (h >= 345):
            return "red"
        if 105 <= h <= 150:
            return "green"
        return "other"

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            C1, h1 = lch(lab[n1]); C2, h2 = lch(lab[n2])
            z1, z2 = zone(h1), zone(h2)
            if {z1, z2} != {"red", "green"}:
                continue
            if C1 <= 40 or C2 <= 40:
                continue
            de = delta_e(lab_de[n1], lab_de[n2])
            out.append(dict(
                color1=n1, hex1=colors[n1], hue1=round(float(h1), 1), C1=round(float(C1), 1),
                color2=n2, hex2=colors[n2], hue2=round(float(h2), 1), C2=round(float(C2), 1),
                dL=round(float(abs(lab[n1][0] - lab[n2][0])), 2),
                deltaE_deuteranopia=round(float(de), 2),
                reference_pure_red_green_deltaE=round(float(ref_de), 2),
                verdict="PURE_RED_GREEN" if de < ref_de else "OK_HUE_SHIFTED",
            ))
    return out
