#!/usr/bin/env python3
"""Figure for the unsupervised report: the two label-free splits in the DiGAS PCA.
A: DiGAS (exon_5k, all ADNI1) PC1 x PC2 coloured by the chr15 ALDH1A2 haplotype group (digas_pca_groups.py)
B: DiGAS without the chr15 block, PC1 x PC2 coloured by the genome-wide (ancestry-like) group (pc1_small_group.csv)
"""
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
B = ROOT / "results/unsupervised_tools/exon_5k"
OUT = ROOT / "results/unsupervised_tools/report"
MUTED, HI, INK, INK2, GRID = "#b5b4ad", "#6250d6", "#0b0b0b", "#52514e", "#d0cfca"

a = pd.read_csv(B / "pc1_groups/subjects.csv")
b = pd.read_csv(B / "no_chr15_block/pc1_small_group.csv")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
panels = [(ax[0], a, a.Group == "left", "A  DiGAS selection: PC1 split",
           "chr15 ALDH1A2 haplotype, no A1 copy", "carries A1"),
          (ax[1], b, b.small_group, "B  DiGAS without the chr15 block: PC1 split",
           "genome-wide (ancestry-like) group", "rest of cohort")]
for axi, d, m, title, lab_hi, lab_lo in panels:
    axi.scatter(d.PC1[~m], d.PC2[~m], s=12, c=MUTED, linewidths=.4, edgecolors="#fcfcfb",
                label=f"{lab_lo} (n={(~m).sum()})")
    axi.scatter(d.PC1[m], d.PC2[m], s=14, c=HI, marker="D", linewidths=.4, edgecolors="#fcfcfb",
                label=f"{lab_hi} (n={m.sum()})")
    axi.set_title(title, loc="left", fontsize=11, color=INK)
    axi.set_xlabel("PC1", color=INK2); axi.set_ylabel("PC2", color=INK2)
    axi.tick_params(colors=INK2, labelsize=8)
    for s in axi.spines.values():
        s.set_color(GRID)
    axi.legend(loc="lower left", fontsize=8, frameon=False, labelcolor=INK)
fig.patch.set_facecolor("#fcfcfb")
fig.tight_layout(); fig.savefig(OUT / "fig_pc1_splits.png", dpi=160)
