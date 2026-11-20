#!/bin/bash
# Result 5: data quality, group differences and statistical power of ADNI1 (no classifier, no permuted labels).
set -e
cd "$(dirname "$0")/.."
.venv/bin/python scripts/analysis/snp_variability.py
