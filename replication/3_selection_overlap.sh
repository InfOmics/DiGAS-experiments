#!/bin/bash
# Result 4: agreement between the regions selected by DiGAS and by SKAT (needs scripts 1 and 2 first).
set -e
cd "$(dirname "$0")/.."
bash scripts/analysis/download_annotation.sh       # public annotation files -> annotation/
.venv/bin/python scripts/analysis/selection_overlap.py
.venv/bin/python scripts/analysis/selection_overlap_annot.py
.venv/bin/python scripts/analysis/selection_overlap_figs.py
