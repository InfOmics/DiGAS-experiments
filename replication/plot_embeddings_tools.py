#!/usr/bin/env python3
"""Embedding figures (PCA, t-SNE, UMAP) of the whole-cohort DiGAS and SKAT selections,
redrawn from the coordinates saved by scripts/analysis/unsupervised_tools.py. Same style as the report figures,
with only the DiGAS and SKAT rows (no reference rows).

  python3 replication/plot_embeddings_tools.py
Outputs: results/unsupervised_tools/report/fig_embeddings_tools_<subset>_<dataset>.png
"""
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "unsupervised_tools"
OUT = BASE / "report"
DIRS = {"ADNI1": "exon_5k", "ADNI2": "exon_5k_ADNI2", "ADNI3": "exon_5k_ADNI3"}
SETS = {"DiGAS": "DiGAS", "SKAT": "SKAT"}
SPACES = {"PCA2D": "PCA", "tSNE2D": "t-SNE", "UMAP2D": "UMAP"}
STATES = ["CN", "MCI", "AD"]
COL = {"CN": "#2a78d6", "MCI": "#eb6834", "AD": "#1baf7a"}
MRK = {"CN": "o", "MCI": "s", "AD": "^"}
SUBSETS = {"CN_AD": "CN / AD", "all": "CN / MCI / AD"}

for ds, d in DIRS.items():
    for sub, sub_title in SUBSETS.items():
        fig, ax = plt.subplots(len(SETS), len(SPACES), figsize=(11, 3.4 * len(SETS)), squeeze=False)
        for i, (s, s_label) in enumerate(SETS.items()):
            for j, (sp, sp_label) in enumerate(SPACES.items()):
                a = ax[i, j]
                C = pd.read_csv(BASE / d / f"coords_{sub}_{s}_{sp}.csv")
                for v in STATES:
                    m = C.State == v
                    if m.any():
                        a.scatter(C.loc[m, f"{sp}_1"], C.loc[m, f"{sp}_2"], s=11, c=COL[v], marker=MRK[v],
                                  alpha=.75, linewidths=.4, edgecolors="#fcfcfb", label=f"{v} (n={m.sum()})")
                a.set_xticks([]); a.set_yticks([])
                for spine in a.spines.values():
                    spine.set_color("#d0cfca")
                if i == 0:
                    a.set_title(sp_label, fontsize=11, color="#0b0b0b")
                if j == 0:
                    a.set_ylabel(s_label, fontsize=11, color="#0b0b0b")
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.suptitle(f"{ds}, exon_5k, whole-cohort selection: {sub_title}", fontsize=12, color="#0b0b0b", y=0.995)
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=len(l), fontsize=9, frameon=False,
                   markerscale=1.6)
        fig.patch.set_facecolor("#fcfcfb")
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        path = OUT / f"fig_embeddings_tools_{sub}_{ds}.png"
        fig.savefig(path, dpi=150); plt.close(fig)
        print("wrote", path.relative_to(ROOT))
