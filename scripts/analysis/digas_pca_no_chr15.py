#!/usr/bin/env python3
"""DiGAS whole-cohort selection (exon_5k) re-analysed without the chr15 ALDH1A2 haplotype block that drives PC1
(see digas_pca_groups.py). Same unsupervised steps as unsupervised_tools.py; DiGAS (full) and Random@DiGAS are
taken from that run's checkpoints for side-by-side comparison.

  .venv/bin/python scripts/analysis/digas_pca_no_chr15.py
Outputs: results/unsupervised_tools/exon_5k/no_chr15_block/
"""
import importlib.util, pickle
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("U", ROOT / "scripts/analysis/unsupervised_tools.py")
U = importlib.util.module_from_spec(spec); spec.loader.exec_module(U)
REGION = "exon_5k"
BASE = ROOT / "results/unsupervised_tools" / REGION
OUT = BASE / "no_chr15_block"; OUT.mkdir(parents=True, exist_ok=True)
DROP = {"ENSG00000128918.15_7", "ENSG00000274719.1_7", "ENSG00000259285.1_6"}   # ALDH1A2 block, chr15:58.2-58.8 Mb
NAME = "DiGAS-chr15"


def main():
    subj, bim, regions = U.C.load_data(ROOT, "ADNI1_QC", ROOT / U.C.DATASETS["ADNI1_QC"], ROOT / "region" / f"{REGION}.txt")
    fam = pd.read_csv(ROOT / "ADNI1_QC.fam", sep=r"\s+", header=None, usecols=[1], names=["IID"], dtype=str)
    lab = subj.set_index(subj.IID.astype(str)).OriginalGroup
    rows = np.where(fam.IID.isin(lab.index))[0]
    y_all = lab.loc[fam.IID.iloc[rows]].values
    bim = bim.reset_index(drop=True)
    r2s = U.region_snps(regions, bim)
    genes = sorted(U.selected_regions(U.SELECTIONS[("DiGAS", REGION)]) & set(r2s))
    kept = [g for g in genes if g not in DROP]
    idx = sorted({i for g in kept for i in r2s[g]})
    subsets = {"all": np.ones(len(rows), bool), "CN_AD": y_all != "MCI"}
    pcs = U.pca_scores(ROOT / "ADNI1_QC", idx, len(fam), rows, subsets)
    print(f"{NAME}: {len(kept)} regions (dropped {len(genes) - len(kept)}), {len(idx)} SNPs")

    # PC1 chromosome shares (loadings via Z: X^T U / s, chunk-free since 5.7k SNPs)
    G = U.read_bed(ROOT / "ADNI1_QC", idx, len(fam))[rows]
    med = np.nanmedian(G, axis=0); ok = ~np.isnan(med)
    X = np.where(np.isnan(G[:, ok]), med[ok], G[:, ok]); X[X == 2] = 1
    keep = X.std(0) > 0; X = X[:, keep]; snp = np.array(idx)[ok][keep]
    Xc = X - X.mean(0)
    _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    evr = s ** 2 / (s ** 2).sum()
    drv = []
    for p in range(3):
        sh = pd.Series(Vt[p] ** 2).groupby(bim.loc[snp, "chr"].values).sum().sort_values(ascending=False)
        drv.append(f"PC{p + 1} ({evr[p]:.3%}): top chr {sh.index[0]} share {sh.iloc[0]:.3f}, "
                   f"2nd chr {sh.index[1]} {sh.iloc[1]:.3f}, max single-SNP share {np.max(Vt[p] ** 2):.4f}")
    print("\n".join(drv))

    rep, clus, embs = drv + [""], [], []
    for sub, m in subsets.items():
        st = {}
        Z, ev, ns = pcs[sub]
        c, e, _ = U.analyse(NAME, 0, Z, ev, ns, y_all[m], 1, st)
        stores = {}
        ref = {}
        for s_name, f in (("DiGAS", "DiGAS_0.pkl"), ("Random@DiGAS", "Random_at_DiGAS_0.pkl")):
            (cc, ee, _), stt = pickle.load(open(BASE / "checkpoints" / f, "rb"))[sub]
            stores.update(stt); ref[s_name] = (cc, ee)
        stores.update(st)
        U.plot(stores, y_all[m], f"ADNI1 {REGION}: DiGAS with / without the chr15 ALDH1A2 block "
               f"({'CN / MCI / AD' if sub == 'all' else 'CN / AD'}; Random@DiGAS draw 1)", OUT / f"embeddings_no_chr15_{sub}.png")
        for (s_name, (cc, ee)) in list(ref.items()) + [(NAME, (c, e))]:
            clus += [dict(r, Subset=sub) for r in cc]; embs += [dict(r, Subset=sub) for r in ee]
        for sp in ("PCA2D", "PCA3D", "tSNE2D", "tSNE3D", "UMAP2D", "UMAP3D"):
            E = st[(NAME, sp)]
            pd.DataFrame(E, columns=[f"{sp}_{i + 1}" for i in range(E.shape[1])]) \
                .assign(IID=fam.IID.iloc[rows][m].values, State=y_all[m], Set=NAME) \
                .to_csv(OUT / f"coords_{sub}_{NAME}_{sp}.csv", index=False)
    clus, embs = pd.DataFrame(clus), pd.DataFrame(embs)
    clus.to_csv(OUT / "clustering.csv", index=False); embs.to_csv(OUT / "embedding.csv", index=False)
    order = ["DiGAS", NAME, "Random@DiGAS"]
    for sub, t, k in (("all", "CN/MCI/AD", 3), ("CN_AD", "CN/AD", 2)):
        a = clus[(clus.Subset == sub) & (clus.Target == t) & (clus.k == k)].pivot_table(
            index="Set", columns="Algorithm", values="ARI").reindex(order).round(4)
        b = embs[(embs.Subset == sub) & (embs.Target == t)].pivot_table(
            index="Set", columns="Space", values="kNN10_excess").reindex(order)[["PCA50", "PCA2D", "tSNE2D", "UMAP2D"]].round(4)
        sil = embs[(embs.Subset == sub) & (embs.Target == t)].pivot_table(
            index="Set", columns="Space", values="Silhouette").reindex(order)[["PCA50", "PCA2D"]].round(4)
        rep += [f"== {sub} ({t}), clustering ARI k={k}", a.to_string(), "-- 10-NN agreement excess", b.to_string(),
                "-- silhouette of true states", sil.to_string(), ""]
    open(OUT / "report.txt", "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
