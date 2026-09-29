# Code Description: Comparative Genomics of Ginsenoside-Relevant Glycoside Hydrolases in Gut Bacteroidota

## Overview
This repository contains the custom analysis scripts, computational pipelines, and configuration files used in the study: *"Comparative genomics of ginsenoside-relevant glycoside hydrolases in gut Bacteroidota: an order-level GH1 detection boundary in Bacteroidales, dual-metric PUL co-localization of GH3, and structure-anchored candidate enzymes for experimental validation."* 

The codebase is designed to ensure full reproducibility of the genomic annotation, statistical evaluation, structural pre-screening, and metatranscriptomic analyses described in the manuscript. No gene re-prediction was performed; all analyses utilized official NCBI RefSeq protein releases.

## Computational Environment & Dependencies
All analyses were executed using the following core software environments and libraries. Exact version numbers and random seeds are documented in the respective script headers and log files.
*   **Programming Languages:** Python 3.11, R 4.2.2
*   **Python Libraries:** NumPy, SciPy, pandas, matplotlib, pyhmmer, Biopython
*   **R Packages:** ape, nlme (used for phylogenetic cross-validation)
*   **External Bioinformatics Tools:** 
    *   HMMER (`hmmsearch`) via `pyhmmer`
    *   MAFFT 7.505 (Multiple sequence alignment, utilizing L-INS-i and FFT-NS-2 strategies)
    *   MUSCLE
    *   DIAMOND (for metatranscriptome blastp searches)
    *   AlphaFold2 & ESMFold (for catalytic geometry modeling)

## Pipeline Architecture
The analytical workflow is modularized into the following key components:

1.  **CAZyme Annotation Module:** Executes `dbCAN.hmm V5-2` scanning against 875 HMM profiles. Implements strict filtering logic (i-evalue < 1e-15, coverage > 0.35, per-target cross-family overlap filtering). Database integrity is verified via MD5 checksums prior to execution.
2.  **PUL Identification & Null Models:** Identifies SusC/SusD pairs using InterPro-released Pfam HMMs. Computes enrichment ratios using three permutation null models (Uniform-window, Model P, and the cluster-preserving circular-shift Model C) with 2,000 permutations each. Fixed random seeds are recorded for exact reproducibility.
3.  **Phylogenetic & Statistical Analysis:** Computes strain-level Kruskal-Wallis tests, Dunn post-hoc tests (Benjamini-Hochberg corrected), Cliff's delta, Fritz-Purvis D, Pagel's lambda, and Blomberg's K.
4.  **Structure-Anchored Candidate Screening:** Applies the pre-specified S1-S4 screening gates. Maps catalytic residues based on structurally fixed columns (e.g., PDB 5Z9S for GH3) rather than propagated database annotations.
5.  **Visualization Engine:** Generates all figures using the Okabe-Ito color-blind-safe palette. Color contrasts are validated under simulated protanopia, deuteranopia, and tritanopia (pairwise ΔE*ab ≥ 16). Every figure output is programmatically linked to a source-data table with a provenance column.

## Reproducibility Notes
*   **Frozen Protocol:** The dbCAN database version (V5-2) and all filtering thresholds were frozen prior to execution.
*   **Audit Trail:** Command logs, software versions, database MD5 checksums, and pre-frozen criteria files are included in the `/logs` and `/config` directories.
*   **Data Provenance:** A number-provenance ledger (Table S8 in the manuscript) maps every quantitative claim to its specific source file and cell within this repository.
