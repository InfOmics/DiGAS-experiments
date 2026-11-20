#!/usr/bin/env python3
"""Figures for the DiGAS/SKAT common-selection report (reads results/selection_overlap/*.csv).
A: share of DiGAS-selected regions also selected by SKAT, by region size (exon_5k, ADNI1/2/3)
B: within-cohort overlap (DiGAS vs SKAT) vs cross-cohort overlap (same tool, two cohorts), observed / chance
"""
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "results/selection_overlap"; OUT = IN / "report"; OUT.mkdir(parents=True, exist_ok=True)
C3 = {"ADNI1": "#2a78d6", "ADNI2": "#eb6834", "ADNI3": "#1baf7a"}     # categorical slots 1-3 (validated order)
INK, INK2, GRID, MUTED, BG = "#0b0b0b", "#52514e", "#d0cfca", "#b5b4ad", "#fcfcfb"


def style(ax):
    ax.tick_params(colors=INK2, labelsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.yaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)


ag = pd.read_csv(IN / "agreement_by_size.csv"); ag = ag[ag.Region == "exon_5k"]
ov = pd.read_csv(IN / "overlap.csv"); xc = pd.read_csv(IN / "cross_cohort.csv")
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw=dict(width_ratios=[1.1, 1]))

# A
bins = ["1", "2", "3-5", "6-10", ">10"]; w = .26
for i, ds in enumerate(C3):
    d = ag[ag.Dataset == ds].set_index("SNPs").reindex(bins)
    x = np.arange(len(bins)) + (i - 1) * w
    ax[0].bar(x, d.Common_of_DiGAS * 100, w - .03, color=C3[ds], label=ds, edgecolor=BG, linewidth=1)
ax[0].set_xticks(range(len(bins))); ax[0].set_xticklabels(bins)
ax[0].set_xlabel("SNPs in the region", color=INK2); ax[0].set_ylabel("% of DiGAS selections also selected by SKAT", color=INK2)
ax[0].set_ylim(0, 100); ax[0].legend(frameon=False, fontsize=8.5, labelcolor=INK)
ax[0].set_title("A  Agreement falls with region size", loc="left", fontsize=11, color=INK)
style(ax[0])

# B: fold over chance (log scale), within-cohort union vs cross-cohort
rows = []
for ds in C3:
    r = ov[(ov.Dataset == ds) & (ov.Region == "exon_5k") & (ov.Contrast == "union")].iloc[0]
    rows.append((ds, r.Common / r.Expected, "within"))
for t in ("DiGAS", "SKAT"):
    for _, r in xc[xc.Set == t].iterrows():
        rows.append((f"{t}\n{r.Pair.replace('ADNI', '').replace('-', '–')}", r.Common / r.Expected, "across"))
lab, fold, kind = zip(*rows)
x = np.arange(len(rows))
cols = [INK2 if k == "within" else MUTED for k in kind]
ax[1].bar(x, fold, .62, color=cols, edgecolor=BG, linewidth=1)
ax[1].axhline(1, color=INK, lw=1, ls=(0, (3, 2)))
for xi, f in zip(x, fold):
    ax[1].text(xi, f * 1.06, f"{f:.1f}×", ha="center", va="bottom", fontsize=7.5, color=INK)
ax[1].set_yscale("log"); ax[1].set_ylim(.4, 12)
ax[1].set_yticks([.5, 1, 2, 5, 10]); ax[1].set_yticklabels(["0.5", "1", "2", "5", "10"]); ax[1].minorticks_off()
ax[1].set_xticks(x); ax[1].set_xticklabels(lab, fontsize=7.5)
ax[1].set_ylabel("observed / chance overlap (log; 1 = chance)", color=INK2)
ax[1].set_title("B  Agreement within vs across cohorts", loc="left", fontsize=11, color=INK)
ax[1].text(1, -.16, "DiGAS vs SKAT,\nsame cohort", ha="center", va="top", fontsize=8, color=INK2, transform=ax[1].get_xaxis_transform())
ax[1].text(6, -.16, "same tool, two cohorts\n(1 = ADNI1, 2 = ADNI2, 3 = ADNI3)", ha="center", va="top", fontsize=8, color=INK2, transform=ax[1].get_xaxis_transform())
style(ax[1])
fig.patch.set_facecolor(BG)
fig.tight_layout(); fig.savefig(OUT / "fig_overlap.png", dpi=160)
print(OUT / "fig_overlap.png")
