#!/bin/bash
# Result 4: DiGAS fold change and rank of the APOE-locus genes (APOE, TOMM40, APOC1) in the whole cohort and in each
# nested training fold. Needs 1_supervised_comparison.sh first (nested DiGAS exon_5k run).
set -e
cd "$(dirname "$0")/.."
.venv/bin/python scripts/analysis/apoe_locus_fc.py
