#!/bin/bash
# Result 3: unsupervised comparison (PCA, clustering, 10-NN agreement) of the whole-cohort exon_5k selections.
# Needs 1_supervised_comparison.sh first (ADNI1 whole-cohort selections).
set -e
cd "$(dirname "$0")/.."
.venv/bin/python scripts/analysis/unsupervised_tools.py --region exon_5k --dataset ADNI1
# ADNI2 and ADNI3: whole-cohort DiGAS and SKAT selections first, then the same analysis
.venv/bin/python scripts/analysis/adni23_selections.py
for ds in ADNI2 ADNI3; do
  .venv/bin/python scripts/analysis/unsupervised_tools.py --region exon_5k --dataset $ds
done
# Figures of the document (DiGAS and SKAT rows only)
.venv/bin/python replication/plot_embeddings_tools.py
# Additional analyses of the full report results/unsupervised_tools/report/unsupervised_DiGAS_vs_SKAT.pdf
.venv/bin/python scripts/analysis/digas_pca_groups.py          # ADNI1: what drives DiGAS PC1
.venv/bin/python scripts/analysis/digas_pca_no_chr15.py        # ADNI1: DiGAS without the chr15 block
.venv/bin/python scripts/analysis/cnad_split_drivers.py --dataset ADNI2   # ADNI2: what drives the CN/AD split
.venv/bin/python scripts/analysis/cnad_split_drivers.py --dataset ADNI1
.venv/bin/python scripts/analysis/unsupervised_report_figs.py
