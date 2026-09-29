# 05_脚本_Scripts · 说明

**本目录内容**：全部绘图与分析脚本（原论文 8 图 + 审稿补数据 1A/1B/1C 分析脚本 + R1–R3 重绘脚本）。

## 子目录

| 子目录 | 内容 |
|---|---|
| `figures_main/` | 原论文绘图轮脚本：plot_fig1–5、plot_suppfig1–3（**suppfig2 为环形新版**）、run_all.py、shared_style.py（样式与路径单一真源）、cvd.py（色盲模拟验证）、make_palette_sheet.py、requirements.txt；`input/outputs_round8/stats/bacteria_midpoint.nwk` 为 SF2 环形重绘 PLOT 阶段所需树文件 |
| `review_1A_结构预筛/` | 1A 分析脚本：geom_lib.py（催化几何测量库）、plot_fig_1A.py（原始绘版）、shared_style.py |
| `review_1B_保簇零模型/` | 1B 分析脚本：run_scans45.py + scan_one.py（45 基因组 dbCAN 重扫）、parse_gff45.py、dbcan_postprocess.py、nullmodels_45.py（model P 复现 + model C 补跑，种子 20260920/20260921）、plot_fig_1B.py（原始绘版）、shared_style.py |
| `review_1C_宏转录组/` | 1C 分析脚本：run_1C_diamond.sh（DIAMOND blastp，--very-sensitive, e<1e-5, id≥60%）、plot_fig_1C.py（原始绘版）、shared_style.py |
| `review_redraw_重绘/` | **R1–R3 期刊级重绘脚本**：review_style.py（统一样式：Arial/Liberation Sans、Okabe–Ito、三格式输出）、build_source_data.py（BUILD 阶段：从上游结果表装配 source_data）、plot_figR1.py / plot_figR2.py / plot_figR3.py（PLOT 阶段）、source_data/（装配好的绘图原始数据） |

## 复现方法（已实测通过）

**原论文 8 图（含环形新版 SF2）**：
```bash
cd figures_main
FIG_INPUT="$PWD/input" python3 run_all.py --plot-only
```
- `--plot-only` 只读 `../main/` 与 `../supplementary/` 下的 source_data.tsv（绘图原始数据，已随包放置）重绘全部 8 图到同目录；实测 8/8 通过。
- 不带 `--plot-only` 的 BUILD 阶段需要原项目 round4–round9 上游目录（本包数据在 `04_数据_Data/`，路径与当时项目布局不同，BUILD 不可直接运行；source_data 已交付，无需重跑）。

**审稿图 R1–R3 重绘**：
```bash
cd review_redraw_重绘
python3 build_source_data.py --src-1a <1A结果表目录> --src-1b <1B结果表目录> --src-1c <1C结果表目录> --out source_data
python3 plot_figR1.py --source-data source_data/figR1_structure_geometry_source_data.tsv --out out
python3 plot_figR2.py --forest source_data/figR2_modelC_forest_source_data.tsv --perlayer source_data/figR2_modelC_perlayer_source_data.tsv --distributions source_data/figR2_modelC_permutation_distributions_source_data.tsv --out out
python3 plot_figR3.py --source-data source_data/figR3_metatranscriptome_expression_source_data.tsv --out out
```
- 1A/1B/1C 结果表目录即 `04_数据_Data/review_1*/结果表/`。
- 输出 SVG+PDF+PNG(300dpi)，数值全部从 source_data 读回、零硬编码。

## 环境

Python 3.11 + NumPy/SciPy/pandas/matplotlib（重绘另需 scipy.stats）；1B 重扫需 pyhmmer、1A 需 Biopython、MAFFT 7.505、1C 需 DIAMOND 2.1.9（版本见 `04_数据_Data/logs/`）。
