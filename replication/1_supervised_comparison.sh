#!/bin/bash
# Results 1 and 2: DiGAS vs SKAT classification, whole-cohort (leaky) and nested selection, ADNI1.
# Run from the repository root. Restart-safe: finished jobs leave a .done file and are skipped.
set -e
cd "$(dirname "$0")/.."
# For ADNI1 exon_5k the queue reuses whole-cohort selections from two fixed paths. Create them if missing.
D=results/fold_pair_scatter_leaky/digas_all_subjects/digas_all
S=results/skat_leaky/skat_all_exon_5k
mkdir -p "$(dirname $D)" "$(dirname $S)"
[ -e "${D}_significant.txt" ] || .venv/bin/python scripts/nested/digas3-linux_validation.py --bfile ADNI1_QC \
    --region region/exon_5k.txt --pheno pheno/adni1pheno.txt --iteration 1000 --out "$D"
[ -e "${S}_significant.txt" ] || .venv/bin/python scripts/analysis/skat_tool.py --bfile ADNI1_QC \
    --region region/exon_5k.txt --pheno pheno/adni1pheno.txt --out "$S"
.venv/bin/python scripts/analysis/tool_comparison_queue.py --datasets ADNI1 \
    --regions exon_5k exon exon_protein_coding --slots 2
# Summary tables: results/tool_comparison/summary_all.csv, summary_folds.csv, summary_paired_nested.csv
.venv/bin/python scripts/analysis/compare_tools.py
