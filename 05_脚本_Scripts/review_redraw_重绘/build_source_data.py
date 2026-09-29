#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_source_data.py — 从重绘上游只读结果表生成三张审稿图的 source_data
=====================================================================
两阶段流程的 BUILD 阶段：本脚本读取补数据上游分析结果表（只读），
生成 Fig R1/R2/R3 绘图阶段直接读取的 *_source_data.tsv；
绘图脚本不再触碰上游文件，图上每个数值都可追溯到 source_data 行。

输入（只读）：
  1A  t1a_catalytic_geometry_screen.tsv + t1a_reference_calibration.tsv
  1B  t3_null45_modelPC_aggregate.tsv / _perlayer.tsv + t3_null45_modelC_distributions.npz
  1C  t1c_metatx_expression.tsv
输出（本目录 source_data/）：
  figR1_structure_geometry_source_data.tsv
  figR2_modelC_forest_source_data.tsv
  figR2_modelC_perlayer_source_data.tsv
  figR2_modelC_permutation_distributions_source_data.tsv
  figR3_metatranscriptome_expression_source_data.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_r1(src_1a: Path, outdir: Path) -> Path:
    scr = pd.read_csv(src_1a / "t1a_catalytic_geometry_screen.tsv", sep="\t")
    ref = pd.read_csv(src_1a / "t1a_reference_calibration.tsv", sep="\t")
    cand = pd.DataFrame({
        "block": "cand",
        "entity": scr["protein_id"],
        "family": scr["family"],
        "rank": scr["rank"],
        "tier": scr["tier"],
        "species": scr["species"],
        "struct_source": scr["struct_source"],
        "pos_acid_base": scr["pos_acid_base"],
        "pos_nucleophile": scr["pos_nucleophile"],
        "cdcd_A": scr["cdcd"],
        "minOO_A": scr["minOO"],
        "caca_A": scr["caca"],
        "plddt_acid_base": scr["plddt_ab"],
        "plddt_nucleophile": scr["plddt_nu"],
        "ref_cdcd_A": scr["ref_cdcd"],
        "ref_minOO_A": scr["ref_minOO"],
        "consistency_label": scr["catalytic_geometry_consistency"],
        "provenance": "t1a_catalytic_geometry_screen.tsv",
    })
    refs = pd.DataFrame({
        "block": "ref",
        "entity": ref["reference"],
        "family": ref["family"],
        "rank": np.nan, "tier": "experimental reference",
        "species": "", "struct_source": "PDB (experimental)",
        "pos_acid_base": ref["pos_acid_base"],
        "pos_nucleophile": ref["pos_nucleophile"],
        "cdcd_A": ref["cdcd"], "minOO_A": ref["minOO"], "caca_A": ref["caca"],
        "plddt_acid_base": np.nan, "plddt_nucleophile": np.nan,
        "ref_cdcd_A": ref["cdcd"], "ref_minOO_A": ref["minOO"],
        "consistency_label": "experimental reference",
        "provenance": "t1a_reference_calibration.tsv",
    })
    out = pd.concat([cand, refs], ignore_index=True)
    p = outdir / "figR1_structure_geometry_source_data.tsv"
    out.to_csv(p, sep="\t", index=False)
    return p


def build_r2(src_1b: Path, outdir: Path) -> list[Path]:
    agg = pd.read_csv(src_1b / "t3_null45_modelPC_aggregate.tsv", sep="\t")
    forest = agg.rename(columns={
        "enrichment_ratio": "ER", "empirical_p_one_sided_ge": "empirical_p"})[
        ["model", "metric", "observed", "n_pairs", "null_mean", "null_sd",
         "null_p02_5", "null_p97_5", "ER", "empirical_p", "nperm", "seed"]]
    forest["provenance"] = "t3_null45_modelPC_aggregate.tsv"
    p1 = outdir / "figR2_modelC_forest_source_data.tsv"
    forest.to_csv(p1, sep="\t", index=False)

    lay = pd.read_csv(src_1b / "t3_null45_modelPC_perlayer.tsv", sep="\t")
    lay = lay.rename(columns={"enrichment_ratio": "ER",
                              "empirical_p_one_sided_ge": "empirical_p"})[
        ["model", "scope", "metric", "observed", "n_pairs", "null_mean",
         "null_sd", "ER", "empirical_p", "nperm", "seed"]]
    lay["provenance"] = "t3_null45_modelPC_perlayer.tsv"
    p2 = outdir / "figR2_modelC_perlayer_source_data.tsv"
    lay.to_csv(p2, sep="\t", index=False)

    z = np.load(src_1b / "t3_null45_modelC_distributions.npz")
    rows = []
    for key in z.files:
        model, metric = key.split("|")[0], key.split("|")[1]
        scope = key.split("|")[2] if key.count("|") == 2 else "all45"
        for i, v in enumerate(z[key]):
            rows.append((model, scope, metric, i, int(v)))
    dist = pd.DataFrame(rows, columns=["model", "scope", "metric", "draw", "null_count"])
    dist["provenance"] = "t3_null45_modelC_distributions.npz"
    p3 = outdir / "figR2_modelC_permutation_distributions_source_data.tsv"
    dist.to_csv(p3, sep="\t", index=False)
    return [p1, p2, p3]


def build_r3(src_1c: Path, outdir: Path) -> Path:
    t = pd.read_csv(src_1c / "t1c_metatx_expression.tsv", sep="\t")
    out = t[["protein_id", "family", "rank", "tier", "species", "pul_window",
             "s4b_flag", "best_pident", "best_qcov", "n_hits",
             "n_ge90_len50", "n_ge95_len50", "metatx_expression_call"]].copy()
    out["provenance"] = "t1c_metatx_expression.tsv"
    p = outdir / "figR3_metatranscriptome_expression_source_data.tsv"
    out.to_csv(p, sep="\t", index=False)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src-1a", required=True, type=Path, help="1A 结果表目录")
    ap.add_argument("--src-1b", required=True, type=Path, help="1B 结果表目录")
    ap.add_argument("--src-1c", required=True, type=Path, help="1C 结果表目录")
    ap.add_argument("--out", required=True, type=Path, help="source_data 输出目录")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    paths = [build_r1(a.src_1a, a.out)]
    paths += build_r2(a.src_1b, a.out)
    paths.append(build_r3(a.src_1c, a.out))
    for p in paths:
        print("built:", p)


if __name__ == "__main__":
    main()
