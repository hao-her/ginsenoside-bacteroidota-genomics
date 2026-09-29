#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py — 一键复现全套交付物
================================
顺序：8 张图的 BUILD（装配各自 source_data.tsv）+ PLOT（只读 source_data 绘图），
最后跑色盲验证与配色样张（palette/）。

用法：
    python run_all.py              # 全量：重建 source_data + 重绘 + 色盲验证
    python run_all.py --plot-only  # 仅从已有 source_data.tsv 重绘（不读上游表）
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

SCRIPTS = ["plot_fig1.py", "plot_fig2.py", "plot_fig3.py", "plot_fig4.py", "plot_fig5.py",
           "plot_suppfig1.py", "plot_suppfig2.py", "plot_suppfig3.py"]


def main() -> int:
    only = "--plot-only" in sys.argv
    rc = 0
    for s in SCRIPTS:
        args = [PY, str(HERE / s)] + (["--plot-only"] if only else [])
        r = subprocess.run(args, cwd=HERE)
        if r.returncode != 0:
            print(f"[run_all] FAILED: {s}")
            rc = 1
    r = subprocess.run([PY, str(HERE / "make_palette_sheet.py")], cwd=HERE)
    if r.returncode != 0:
        print("[run_all] FAILED: make_palette_sheet.py")
        rc = 1
    print("[run_all] done" if rc == 0 else "[run_all] completed with errors")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
