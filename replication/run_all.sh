#!/bin/bash
# Runs all analyses in order. The supervised step takes several days on a 16 GB machine.
set -e
d="$(dirname "$0")"
bash "$d/1_supervised_comparison.sh"
bash "$d/2_unsupervised_analysis.sh"
bash "$d/3_selection_overlap.sh"
bash "$d/4_apoc1_locus.sh"
bash "$d/5_power_analysis.sh"
