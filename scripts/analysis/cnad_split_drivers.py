#!/usr/bin/env python3
"""What drives the CN-vs-AD clustering of the DiGAS whole-cohort selection (exon_5k, CN+AD subjects)?

In ADNI2 k-means (k=2) on the 50 PCs of the DiGAS regions' SNPs gives ARI 0.70 vs CN/AD (SKAT 0.02); in ADNI1
ARI 0.33 (k=3). Same coding and PCA as unsupervised_tools.py (carrier 0/1, median imputation, CN+AD only).

1. reproduce the clustering; cluster x state / sex / site; per-subject missingness and heterozygosity (batch check)
2. concentration: how PC1 (the PC most tied to state) spreads over SNPs, regions, chromosomes
3. per-SNP CN-vs-AD association and its link to the PC loadings (a label-built axis vs a haplotype block)
4. ablation: drop the top-k regions by PC share, and random halves of the regions -> ARI
5. contrast subsets: regions selected for CN-vs-AD only / CN-vs-MCI only, for DiGAS and SKAT
6. ancestry: correlation of the cluster axis with the top PCs of all genotyped regions (AllRegions)

  .venv/bin/python scripts/analysis/cnad_split_drivers.py --dataset ADNI2
Outputs: results/unsupervised_tools/exon_5k[_<dataset>]/cnad_split/
"""
import argparse, importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu, spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("U", ROOT / "scripts/analysis/unsupervised_tools.py")
U = importlib.util.module_from_spec(spec); spec.loader.exec_module(U)
REGION = "exon_5k"


def carriers_idx(G, idx):
    med = np.nanmedian(G, axis=0); ok = ~np.isnan(med)
    X = np.where(np.isnan(G[:, ok]), med[ok], G[:, ok]); X[X == 2] = 1
    keep = X.std(axis=0) > 0
    return X[:, keep], np.asarray(idx)[ok][keep]


def pcs(X, n=50):
    Xc = X - X.mean(0)
    u, s, vt = np.linalg.svd(Xc, full_matrices=False)
    return u[:, :n] * s[:n], vt[:n], s ** 2 / (s ** 2).sum()


def km_ari(Z, y, ks=(2, 3)):
    return {k: ARI(y, KMeans(k, n_init=20, random_state=1).fit_predict(Z)) for k in ks}


def selected(prefix, contrasts):
    sig = pd.read_csv(f"{prefix}_significant.txt", sep=r"\s+", dtype=str)
    sig = sig[pd.to_numeric(sig.Pvalue, errors="coerce") <= 0.05]
    pair = sig.Cat1 + "-" + sig.Cat2
    want = set()
    for c in contrasts:
        a, b = c.split("-"); want |= {f"{a}-{b}", f"{b}-{a}"}
    return set(sig[pair.isin(want)].Gene.dropna())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dataset", default="ADNI2")
    ds = ap.parse_args().dataset
    out = ROOT / "results/unsupervised_tools" / (REGION if ds == "ADNI1" else f"{REGION}_{ds}") / "cnad_split"
    out.mkdir(parents=True, exist_ok=True)
    pfx = ROOT / f"{ds}_QC"
    subj, bim, regions = U.C.load_data(ROOT, f"{ds}_QC", ROOT / U.C.DATASETS[f"{ds}_QC"], ROOT / "region" / f"{REGION}.txt")
    fam = pd.read_csv(f"{pfx}.fam", sep=r"\s+", header=None, usecols=[1, 4], names=["IID", "Sex"], dtype=str)
    lab = subj.set_index(subj.IID.astype(str)).OriginalGroup
    rows_all = np.where(fam.IID.isin(lab.index))[0]
    st = lab.loc[fam.IID.iloc[rows_all]].values
    rows = rows_all[np.isin(st, ["CN", "AD"])]
    info = fam.iloc[rows].reset_index(drop=True)
    info["State"] = lab.loc[info.IID].values
    info["Sex"] = info.Sex.map({"1": "M", "2": "F"})
    info["Site"] = info.IID.str[:3]
    y = info.State.values; yb = (y == "AD").astype(int)
    bim = bim.reset_index(drop=True)
    r2s = U.region_snps(regions, bim); n_all = len(fam)
    snp2reg = {}
    for r, ii in r2s.items():
        for i in ii:
            snp2reg.setdefault(i, []).append(r)

    def build(genes):
        idx = sorted({i for g in genes for i in r2s[g]})
        X, snp = carriers_idx(U.read_bed(pfx, idx, n_all)[rows], idx)
        return X, snp

    sel = {t: U.selected_regions(U.selection_prefix(t, REGION, ds)) & set(r2s) for t in ("DiGAS", "SKAT")}
    genes = sorted(sel["DiGAS"])
    X, snp = build(genes)
    Z, Vt, evr = pcs(X)
    c2 = KMeans(2, n_init=20, random_state=1).fit_predict(Z)
    c3 = KMeans(3, n_init=20, random_state=1).fit_predict(Z)
    best_k = 2 if ARI(y, c2) >= ARI(y, c3) else 3
    c = c2 if best_k == 2 else c3
    info["Cluster"] = c
    for i in range(3):
        info[f"PC{i + 1}"] = Z[:, i]
    rep = [f"{ds} {REGION}, CN+AD: {len(rows)} subjects ({(y == 'CN').sum()} CN, {(y == 'AD').sum()} AD); "
           f"DiGAS {len(genes)} regions, {X.shape[1]} SNPs",
           f"k-means on 50 PCs: ARI k=2 {ARI(y, c2):.3f}, k=3 {ARI(y, c3):.3f}  (using k={best_k} below)", ""]

    # PC most tied to state
    from scipy.stats import kruskal
    kw = [kruskal(Z[y == "CN", i], Z[y == "AD", i]).pvalue for i in range(Z.shape[1])]
    rep.append("PCs vs state (Mann-Whitney/KW p) and variance explained, top 5 by p:")
    for i in np.argsort(kw)[:5]:
        rep.append(f"  PC{i + 1}: p = {kw[i]:.2g}, var {evr[i]:.3%}")
    pk = int(np.argmin(kw))
    rep += [f"cluster vs sign of PC{pk + 1}: ARI {ARI(c, Z[:, pk] > np.median(Z[:, pk])):.3f}; "
            f"state vs PC{pk + 1}>median: ARI {ARI(y, Z[:, pk] > np.median(Z[:, pk])):.3f}", ""]

    # 1. composition and batch
    for col in ("State", "Sex", "Site"):
        t = pd.crosstab(info[col], info.Cluster); p = chi2_contingency(t)[1]
        t.to_csv(out / f"cluster_vs_{col.lower()}.csv")
        rep += [f"Cluster x {col}: chi2 p = {p:.3g}"] + ([t.to_string()] if col != "Site" else [f"({t.shape[0]} sites)"]) + [""]
    nb = len(pd.read_csv(f"{pfx}.bim", sep=r"\s+", header=None, usecols=[0]))
    miss = np.zeros(len(rows)); het = np.zeros(len(rows)); nobs = np.zeros(len(rows))
    for i in range(0, nb, 20000):
        Gq = U.read_bed(pfx, range(i, min(i + 20000, nb)), n_all)[rows]
        m = np.isnan(Gq); miss += m.sum(1); nobs += (~m).sum(1); het += (Gq == 1).sum(1)
    info["Missing_rate"] = miss / nb; info["Het_rate"] = het / nobs
    for col in ("Missing_rate", "Het_rate"):
        g = [info.loc[info.Cluster == k, col] for k in sorted(info.Cluster.unique())]
        rep.append(f"{col}: medians by cluster " + ", ".join(f"{x.median():.5f}" for x in g) +
                   (f", Mann-Whitney p = {mannwhitneyu(*g).pvalue:.3g}" if len(g) == 2 else ""))
    rep.append("")

    # 2./3. concentration and association
    L = pd.DataFrame({"SNP": bim.iloc[snp, 1].values, "chr": bim.loc[snp, "chr"].values, "bp": bim.loc[snp, "bp"].values,
                      "Regions": [";".join(snp2reg.get(i, [])) for i in snp], "Loading2": Vt[pk] ** 2,
                      "Loading": Vt[pk]})
    fcn, fad = X[yb == 0].mean(0), X[yb == 1].mean(0)
    L["Carrier_CN"], L["Carrier_AD"], L["Diff_AD_CN"] = fcn, fad, fad - fcn
    L["Fisher_p"] = [fisher_exact([[X[yb == 1, j].sum(), yb.sum() - X[yb == 1, j].sum()],
                                   [X[yb == 0, j].sum(), (1 - yb).sum() - X[yb == 0, j].sum()]])[1] for j in range(X.shape[1])]
    L = L.sort_values("Loading2", ascending=False); L.to_csv(out / "snps.csv", index=False)
    cum = L.Loading2.cumsum().values
    eff = 1 / (L.Loading2 ** 2).sum()
    rep += [f"PC{pk + 1} concentration: top-10 SNPs hold {cum[9]:.1%}, top-100 {cum[99]:.1%}; "
            f"SNPs needed for 50% {np.searchsorted(cum, .5) + 1}, 90% {np.searchsorted(cum, .9) + 1} of {len(L)}; "
            f"effective number of SNPs {eff:.0f}"]
    bychr = L.groupby("chr").agg(SNPs=("SNP", "size"), Share=("Loading2", "sum")).sort_values("Share", ascending=False)
    bychr.to_csv(out / "by_chr.csv")
    rep += [f"top chromosome share {bychr.Share.iloc[0]:.1%} (chr {bychr.index[0]}, {bychr.SNPs.iloc[0]} SNPs = "
            f"{bychr.SNPs.iloc[0] / len(L):.1%} of SNPs)"]
    rr = L.assign(Region=L.Regions.str.split(";")).explode("Region").groupby("Region") \
        .agg(chr=("chr", "first"), SNPs=("SNP", "size"), Share=("Loading2", "sum"),
             min_Fisher_p=("Fisher_p", "min")).sort_values("Share", ascending=False)
    rr.to_csv(out / "regions.csv")
    rc = rr.Share.cumsum().values / rr.Share.sum()
    rep += [f"regions needed for 50% of PC{pk + 1}: {np.searchsorted(rc, .5) + 1}, 90%: {np.searchsorted(rc, .9) + 1} "
            f"of {len(rr)}", "top 10 regions:", rr.head(10).round(4).to_string(), ""]
    sgn = np.sign(np.corrcoef(Z[:, pk], yb)[0, 1])
    rho = spearmanr(L.Loading * sgn, L.Diff_AD_CN)
    rep += [f"SNPs with Fisher p < 0.05 for CN vs AD: {(L.Fisher_p < .05).mean():.1%} (null ~5%), p < 1e-3: "
            f"{(L.Fisher_p < 1e-3).sum()}, p < 1e-6: {(L.Fisher_p < 1e-6).sum()}",
            f"Spearman(signed PC{pk + 1} loading, AD-CN carrier difference) = {rho.correlation:.3f} (p {rho.pvalue:.2g})",
            f"label-aligned score: corr(PC{pk + 1}, sum over SNPs of sign(AD-CN diff)*carrier) = "
            f"{abs(np.corrcoef(Z[:, pk], (X - X.mean(0)) @ np.sign(fad - fcn))[0, 1]):.3f}", ""]

    # 4. ablation
    rep.append(f"ablation (k-means ARI vs CN/AD, k=2 / k=3):")
    abl = []
    order = list(rr.index)
    for k in (0, 1, 10, 50, 100, 300):
        g = [r for r in genes if r not in set(order[:k])]
        Xa, _ = build(g); a = km_ari(pcs(Xa)[0], y)
        abl.append(dict(Removed=f"top {k}", Regions=len(g), ARI_k2=a[2], ARI_k3=a[3]))
    rng = np.random.default_rng(1)
    for d in range(5):
        g = list(rng.choice(genes, len(genes) // 2, replace=False))
        Xa, _ = build(g); a = km_ari(pcs(Xa)[0], y)
        abl.append(dict(Removed=f"random half {d + 1}", Regions=len(g), ARI_k2=a[2], ARI_k3=a[3]))

    # 5. contrast subsets for both tools
    for t in ("DiGAS", "SKAT"):
        pf = U.selection_prefix(t, REGION, ds)
        for name, cons in (("CN-AD only", ["CN-AD"]), ("CN-MCI only", ["CN-MCI"]), ("CN-AD + CN-MCI (as used)", ["CN-AD", "CN-MCI"])):
            g = sorted(selected(pf, cons) & set(r2s))
            Xa, _ = build(g); a = km_ari(pcs(Xa)[0], y)
            abl.append(dict(Removed=f"{t} {name}", Regions=len(g), ARI_k2=a[2], ARI_k3=a[3]))
    abl = pd.DataFrame(abl); abl.to_csv(out / "ablation.csv", index=False)
    rep += [abl.round(3).to_string(index=False), ""]

    # 6. ancestry axis
    allidx = sorted({i for r in r2s for i in r2s[r]})
    Za = U.pca_scores(pfx, allidx, n_all, rows, {"cnad": np.ones(len(rows), bool)}, n_comp=10)["cnad"][0]
    rep.append(f"|corr| of DiGAS PC{pk + 1} with AllRegions PC1..PC5: " +
               ", ".join(f"{abs(np.corrcoef(Z[:, pk], Za[:, i])[0, 1]):.3f}" for i in range(5)))
    rep.append("AllRegions PC1..PC5 vs state (Mann-Whitney p): " +
               ", ".join(f"{mannwhitneyu(Za[yb == 0, i], Za[yb == 1, i]).pvalue:.2g}" for i in range(5)))
    a = km_ari(Za, y); rep.append(f"AllRegions k-means ARI: k=2 {a[2]:.3f}, k=3 {a[3]:.3f}")
    info.to_csv(out / "subjects.csv", index=False)
    open(out / "report.txt", "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
