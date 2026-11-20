#!/usr/bin/env python3
"""DiGAS vs SKAT compared without a classifier: unsupervised learning on the SNPs of the regions each tool selects
on the WHOLE ADNI1 cohort, validated against the known disease state (external validation).

Feature sets (all ADNI1 subjects, same genotype coding as the classification pipeline: carrier 0/1 of A1,
median imputation):
  DiGAS          regions selected by DiGAS on all subjects (CN-vs-MCI / CN-vs-AD contrasts, p <= 0.05)
  SKAT           regions selected by SKAT on all subjects (same rule)
  Random@DiGAS   label-free reference: as many random regions of the region file as DiGAS selected (R draws)
  Random@SKAT    same, matched to the SKAT count
  AllRegions     every SNP inside the region file
The references use no labels at all: they show how much state structure appears in genotypes of the same kind
without any supervised selection. Both tools selected with the labels of these very subjects, so part of any
separation is built in by the selection; the comparison is DiGAS vs SKAT, each read against its reference.

Unsupervised steps (per feature set, per subject set: all CN/MCI/AD, and CN+AD only):
  PCA (centred, unscaled) -> 50 PCs -> clustering: k-means, Ward agglomerative, Gaussian mixture, k = 2 and 3
  external indices vs disease state: ARI, AMI, NMI, V-measure, purity, Fowlkes-Mallows
  2D/3D embeddings: PCA, t-SNE, UMAP; silhouette of the true states and 10-NN state agreement (with its
  expectation under no structure, (n_c - 1)/(n - 1), computed analytically)
  Kruskal-Wallis of the top 10 PCs across states

  .venv/bin/python scripts/analysis/unsupervised_tools.py [--region exon_5k] [--draws 10]
Outputs: results/unsupervised_tools/<region>/
"""
import argparse, importlib.util, os, pickle, sys, time
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "2")

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn import metrics as M
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("C", ROOT / "scripts/nested/digas_classification_nested_10fold.py")
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
SELECTIONS = {  # whole-cohort selections (prefix of <prefix>_significant.txt)
    ("DiGAS", "exon_5k"): ROOT / "results/fold_pair_scatter_leaky/digas_all_subjects/digas_all",
    ("SKAT", "exon_5k"): ROOT / "results/skat_leaky/skat_all_exon_5k",
}


def selection_prefix(tool, region, ds="ADNI1"):
    """ADNI1: the existing whole-cohort selections; ADNI2/ADNI3: the queue's selection_all/<tool>_<region>_<ds>."""
    if ds == "ADNI1":
        return SELECTIONS[(tool, region)]
    return ROOT / "results/tool_comparison/selection_all" / f"{tool}_{region}_{ds}"
STATES = ["CN", "MCI", "AD"]
COL = {"CN": "#2a78d6", "MCI": "#eb6834", "AD": "#1baf7a"}          # categorical slots 1-3 (validated order)
MRK = {"CN": "o", "MCI": "s", "AD": "^"}                             # secondary encoding


def mem_free_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) / 1024 ** 2


def read_bed(prefix, snp_idx, n):
    """A1 allele counts (0/1/2, NaN missing) for the given .bim row indices; plink SNP-major .bed."""
    nb = (n + 3) // 4
    lut = np.array([2, np.nan, 1, 0], dtype=np.float32)              # 00 hom A1, 01 miss, 10 het, 11 hom A2
    out = np.empty((n, len(snp_idx)), dtype=np.float32)
    with open(f"{prefix}.bed", "rb") as f:
        assert f.read(3) == b"\x6c\x1b\x01"
        mm = np.memmap(f"{prefix}.bed", dtype=np.uint8, mode="r", offset=3)
        for j, s in enumerate(snp_idx):
            b = mm[s * nb:(s + 1) * nb]
            codes = np.stack([(b >> sh) & 3 for sh in (0, 2, 4, 6)], axis=1).ravel()[:n]
            out[:, j] = lut[codes]
    return out


def region_snps(regions, bim):
    """region -> list of .bim row indices, same rule as the pipeline (chr equal, start < bp < end)."""
    res = {}
    for c, g in regions.groupby("chr", sort=False):
        b = bim[bim.chr == c].sort_values("bp", kind="stable")
        bp, idx = b.bp.values, b.index.values
        lo = np.searchsorted(bp, g.start.values, side="right"); hi = np.searchsorted(bp, g.end.values, side="left")
        for gene, l, h in zip(g.Gene.values, lo, hi):
            if h > l:
                res.setdefault(gene, []).extend(idx[l:h].tolist())
    return res


def selected_regions(prefix):
    sig = pd.read_csv(f"{prefix}_significant.txt", sep=r"\s+", dtype=str)
    sig["Pvalue"] = pd.to_numeric(sig.Pvalue, errors="coerce")
    sig = sig[sig.Pvalue <= 0.05]
    keep = ((sig.Cat1 == "CN") & sig.Cat2.isin(["AD", "MCI"])) | ((sig.Cat2 == "CN") & sig.Cat1.isin(["AD", "MCI"]))
    return set(sig[keep].Gene.dropna().astype(str))


def carriers(G):
    """pipeline coding: median imputation, then 2 -> 1 (carrier of A1); drop constant columns."""
    med = np.nanmedian(G, axis=0)
    ok = ~np.isnan(med)
    G = G[:, ok]; med = med[ok]
    G = np.where(np.isnan(G), med, G); G[G == 2] = 1
    return G[:, G.std(axis=0) > 0]


def pca_scores(prefix, idx, n, rows, masks, n_comp=50, chunk=10000):
    """exact PCA scores per subject subset without holding the full SNP matrix: the centred Gram matrix
    (subjects x subjects) is accumulated over SNP chunks (carrier coding per subset, as carriers()), then
    eigendecomposed. Z = U*S equals PCA(X).transform(X) up to the sign of each component."""
    K = {k: np.zeros((m.sum(), m.sum())) for k, m in masks.items()}
    n_snp = {k: 0 for k in masks}
    for i in range(0, len(idx), chunk):
        G = read_bed(prefix, idx[i:i + chunk], n)[rows]
        for k, m in masks.items():
            X = carriers(G[m]).astype(np.float64)
            X -= X.mean(axis=0)
            K[k] += X @ X.T; n_snp[k] += X.shape[1]
        del G
    out = {}
    for k, Kk in K.items():
        w, U = np.linalg.eigh(Kk)
        w, U = w[::-1], U[:, ::-1]
        nc = min(n_comp, n_snp[k], Kk.shape[0] - 1)
        w = np.clip(w, 0, None)
        out[k] = (U[:, :nc] * np.sqrt(w[:nc]), w[:nc] / w.sum(), n_snp[k])
    return out


def mem_guard(min_gb=2.0):
    if mem_free_gb() < min_gb:
        sys.exit(f"MemAvailable {mem_free_gb():.1f} GB < {min_gb} GB: stopping (checkpoints kept, rerun to resume)")


def purity(y, c):
    t = pd.crosstab(c, y)
    return t.max(axis=1).sum() / t.values.sum()


def knn_agreement(Z, y, k=10):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Z)
    idx = nn.kneighbors(Z, return_distance=False)[:, 1:]
    y = np.asarray(y)
    obs = (y[idx] == y[:, None]).mean()
    cnt = pd.Series(y).value_counts()
    exp = np.mean([(cnt[v] - 1) / (len(y) - 1) for v in y])
    return obs, exp


def analyse(name, draw, Z, evr, n_snps, y3, seed, emb_store):
    """all analyses for one feature set given its PCA scores Z (<= 50 PCs) and explained variance ratios;
    y3 = CN/MCI/AD labels of the rows."""
    y2 = np.where(y3 == "CN", "CN", "AD+MCI")
    targets = {"CN/MCI/AD": y3} if (y3 == "MCI").any() else {"CN/AD": y3}
    if (y3 == "MCI").any():
        targets["CN vs AD+MCI"] = y2
    clu, emb = [], []
    for k in (2, 3):
        algos = {"KMeans": KMeans(k, n_init=20, random_state=seed).fit_predict(Z),
                 "Ward": AgglomerativeClustering(k, linkage="ward").fit_predict(Z),
                 "GMM": GaussianMixture(k, covariance_type="full", n_init=5, random_state=seed).fit(Z).predict(Z)}
        for a, c in algos.items():
            for t, y in targets.items():
                clu.append({"Algorithm": a, "k": k, "Target": t, "ARI": M.adjusted_rand_score(y, c),
                            "AMI": M.adjusted_mutual_info_score(y, c), "NMI": M.normalized_mutual_info_score(y, c),
                            "V_measure": M.v_measure_score(y, c), "Purity": purity(y, c),
                            "FMI": M.fowlkes_mallows_score(y, c),
                            "Cluster_sizes": "/".join(map(str, np.bincount(c)))})
    spaces = {"PCA50": Z}
    for d in (2, 3):
        spaces[f"PCA{d}D"] = Z[:, :d]
        spaces[f"tSNE{d}D"] = TSNE(d, perplexity=30, init="pca", random_state=seed).fit_transform(Z)
        import umap
        spaces[f"UMAP{d}D"] = umap.UMAP(n_components=d, n_neighbors=15, min_dist=0.1,
                                        random_state=seed).fit_transform(Z)
    for s, E in spaces.items():
        for t, y in targets.items():
            obs, exp = knn_agreement(E, y)
            emb.append({"Space": s, "Target": t, "Silhouette": M.silhouette_score(E, y),
                        "kNN10_agreement": obs, "kNN10_expected": exp, "kNN10_excess": obs - exp})
        if draw == 0 and s != "PCA50":
            emb_store[(name, s)] = E
    pcs = [{"PC": i + 1, "Var_explained": evr[i],
            "KW_p_states": kruskal(*[Z[y3 == v, i] for v in np.unique(y3)]).pvalue} for i in range(min(10, Z.shape[1]))]
    tag = lambda rows: [dict(Set=name, Draw=draw, N_subjects=len(y3), N_SNPs=n_snps, **r) for r in rows]
    return tag(clu), tag(emb), tag(pcs)


def plot(emb_store, y3, title, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sets = list(dict.fromkeys(k[0] for k in emb_store)); spaces = ["PCA2D", "tSNE2D", "UMAP2D"]
    fig, ax = plt.subplots(len(sets), 3, figsize=(11, 3.4 * len(sets)), squeeze=False)
    for i, s in enumerate(sets):
        for j, sp in enumerate(spaces):
            a, E = ax[i, j], emb_store[(s, sp)]
            for v in STATES:
                m = y3 == v
                if m.any():
                    a.scatter(E[m, 0], E[m, 1], s=11, c=COL[v], marker=MRK[v], alpha=.75, linewidths=.4,
                              edgecolors="#fcfcfb", label=f"{v} (n={m.sum()})")
            a.set_xticks([]); a.set_yticks([])
            for sp_ in a.spines.values():
                sp_.set_color("#d0cfca")
            if i == 0:
                a.set_title(sp.replace("2D", ""), fontsize=11, color="#0b0b0b")
            if j == 0:
                a.set_ylabel(s, fontsize=11, color="#0b0b0b")
    ax[0, 2].legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=12, color="#0b0b0b"); fig.patch.set_facecolor("#fcfcfb")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", default="exon_5k")
    ap.add_argument("--dataset", default="ADNI1", choices=["ADNI1", "ADNI2", "ADNI3"])
    ap.add_argument("--draws", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if mem_free_gb() < 3:
        sys.exit("less than 3 GB free RAM: not starting (pipeline jobs running)")
    ds = a.dataset; bfile = ROOT / f"{ds}_QC"
    out = ROOT / "results" / "unsupervised_tools" / (a.region if ds == "ADNI1" else f"{a.region}_{ds}"); out.mkdir(parents=True, exist_ok=True)
    subj, bim, regions = C.load_data(ROOT, f"{ds}_QC", ROOT / C.DATASETS[f"{ds}_QC"], ROOT / "region" / f"{a.region}.txt")
    fam = pd.read_csv(f"{bfile}.fam", sep=r"\s+", header=None, usecols=[1], names=["IID"], dtype=str)
    lab = subj.set_index(subj.IID.astype(str)).OriginalGroup
    rows = np.where(fam.IID.isin(lab.index))[0]
    y_all = lab.loc[fam.IID.iloc[rows]].values
    bim = bim.reset_index(drop=True)
    r2s = region_snps(regions, bim)
    t0 = time.time()

    sets = {}
    for tool in ("DiGAS", "SKAT"):
        genes = selected_regions(selection_prefix(tool, a.region, ds)) & set(r2s)
        sets[tool] = [sorted(genes)]
    rng = np.random.default_rng(a.seed)
    pool = sorted(r2s)
    for tool in ("DiGAS", "SKAT"):
        n = len(sets[tool][0])
        sets[f"Random@{tool}"] = [list(rng.choice(pool, n, replace=False)) for _ in range(a.draws)]
    sets["AllRegions"] = [pool]
    pd.DataFrame([{"Set": s, "Draw": d, "N_regions": len(g)} for s, gs in sets.items() for d, g in enumerate(gs)]) \
        .to_csv(out / "feature_sets.csv", index=False)

    n_all = len(fam)
    subsets = {"all": np.ones(len(rows), bool), "CN_AD": y_all != "MCI"}
    res = {k: ([], [], []) for k in subsets}
    stores = {k: {} for k in subsets}
    ck = out / "checkpoints"; ck.mkdir(exist_ok=True)       # one pickle per (set, draw): a rerun resumes
    for s, gs in sets.items():
        for d, genes in enumerate(gs):
            f = ck / f"{s.replace('@', '_at_')}_{d}.pkl"
            if f.exists():
                part = pickle.load(open(f, "rb"))
            else:
                mem_guard()
                idx = sorted({i for g in genes for i in r2s[g]})
                pcs_by_sub = pca_scores(bfile, idx, n_all, rows, subsets)
                part = {}
                for sub, m in subsets.items():
                    st = {}
                    Z, evr, ns = pcs_by_sub[sub]
                    part[sub] = (analyse(s, d, Z, evr, ns, y_all[m], a.seed + d, st), st)
                pickle.dump(part, open(f.with_suffix(".tmp"), "wb")); f.with_suffix(".tmp").rename(f)
                print(f"[{time.time() - t0:6.0f}s] {s} draw {d}: {len(genes)} regions, {len(idx)} SNPs, "
                      f"free RAM {mem_free_gb():.1f} GB", flush=True)
            for sub in subsets:
                (c, e, p), st = part[sub]
                for acc, x in zip(res[sub], (c, e, p)):
                    acc.extend(x)
                stores[sub].update(st)

    for sub in subsets:
        clu, emb, pcs = (pd.DataFrame(x) for x in res[sub])
        clu.to_csv(out / f"clustering_{sub}.csv", index=False)
        emb.to_csv(out / f"embedding_{sub}.csv", index=False)
        pcs.to_csv(out / f"pcs_{sub}.csv", index=False)
        for (s, sp), E in stores[sub].items():
            if sp.endswith("3D") or sp.endswith("2D"):
                pd.DataFrame(E, columns=[f"{sp}_{i + 1}" for i in range(E.shape[1])]) \
                    .assign(IID=fam.IID.iloc[rows][subsets[sub]].values, State=y_all[subsets[sub]], Set=s) \
                    .to_csv(out / f"coords_{sub}_{s}_{sp}.csv", index=False)
        ylab = y_all[subsets[sub]]
        plot(stores[sub], ylab, f"{ds} {a.region}: {'CN / MCI / AD' if sub == 'all' else 'CN / AD'} "
             f"(Random@ = draw 1 of {a.draws})", out / f"embeddings_{sub}.png")

        # summary: tools vs reference mean +- sd over draws
        def summ(df, keys, cols):
            g = df.groupby(["Set"] + keys)[cols]
            return g.mean().round(4).join(g.std().round(4).add_suffix("_sd"))
        summ(clu, ["Algorithm", "k", "Target"], ["ARI", "AMI", "Purity"]).to_csv(out / f"summary_clustering_{sub}.csv")
        summ(emb, ["Space", "Target"], ["Silhouette", "kNN10_excess"]).to_csv(out / f"summary_embedding_{sub}.csv")
    print(f"done in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
