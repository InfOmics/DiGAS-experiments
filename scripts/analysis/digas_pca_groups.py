#!/usr/bin/env python3
"""What separates the two groups seen along PC1 of the DiGAS whole-cohort selection (exon_5k, all ADNI1 subjects)?

1. groups: 2-component GMM on PC1 (same PCA as unsupervised_tools.py)
2. composition vs disease state, sex, site; per-subject missingness and heterozygosity (batch check)
3. which SNPs/regions/chromosomes drive PC1 (squared loadings), per-SNP carrier frequency difference between groups
4. how few SNPs discriminate the groups (cross-validated logistic regression on the top-k SNPs, and on one region)
5. whether the same subjects split along PC1 of the other feature sets (Random@DiGAS draw 0, SKAT, AllRegions)

  .venv/bin/python scripts/analysis/digas_pca_groups.py
Outputs: results/unsupervised_tools/exon_5k/pc1_groups/
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import StratifiedKFold, cross_val_score

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("U", ROOT / "scripts/analysis/unsupervised_tools.py")
U = importlib.util.module_from_spec(spec); spec.loader.exec_module(U)
REGION = "exon_5k"
OUT = ROOT / "results/unsupervised_tools" / REGION / "pc1_groups"; OUT.mkdir(parents=True, exist_ok=True)
PFX = ROOT / "ADNI1_QC"


def split_pc1(Z, seed=1):
    g = GaussianMixture(2, n_init=10, random_state=seed).fit(Z[:, :1]).predict(Z[:, :1])
    if Z[g == 0, 0].mean() > Z[g == 1, 0].mean():   # 0 = left (lower PC1), 1 = right
        g = 1 - g
    return g


def main():
    subj, bim, regions = U.C.load_data(ROOT, "ADNI1_QC", ROOT / U.C.DATASETS["ADNI1_QC"], ROOT / "region" / f"{REGION}.txt")
    fam = pd.read_csv(f"{PFX}.fam", sep=r"\s+", header=None, usecols=[1, 4], names=["IID", "Sex"], dtype=str)
    lab = subj.set_index(subj.IID.astype(str)).OriginalGroup
    rows = np.where(fam.IID.isin(lab.index))[0]
    info = fam.iloc[rows].reset_index(drop=True)
    info["State"] = lab.loc[info.IID].values
    info["Sex"] = info.Sex.map({"1": "M", "2": "F"})
    info["Site"] = info.IID.str[:3]
    bim = bim.reset_index(drop=True)
    r2s = U.region_snps(regions, bim)
    n_all = len(fam)
    one = {"all": np.ones(len(rows), bool)}

    # --- DiGAS PCA with loadings (5.8k SNPs: the full matrix is small)
    genes = sorted(U.selected_regions(U.SELECTIONS[("DiGAS", REGION)]) & set(r2s))
    idx = sorted({i for g in genes for i in r2s[g]})
    G = U.read_bed(PFX, idx, n_all)[rows]
    med = np.nanmedian(G, axis=0); ok = ~np.isnan(med)
    Gi = np.where(np.isnan(G[:, ok]), med[ok], G[:, ok]); Gi[Gi == 2] = 1
    keep = Gi.std(axis=0) > 0
    X = Gi[:, keep]; snp = np.array(idx)[ok][keep]
    Xc = X - X.mean(axis=0)
    Uu, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    Z = Uu[:, :10] * s[:10]
    evr = s ** 2 / (s ** 2).sum()
    g = split_pc1(Z)
    info["PC1"], info["PC2"], info["Group"] = Z[:, 0], Z[:, 1], np.where(g == 1, "right", "left")
    info.to_csv(OUT / "subjects.csv", index=False)

    rep = [f"DiGAS {REGION}: {len(genes)} regions, {X.shape[1]} SNPs, {len(rows)} subjects",
           f"PC1 explains {evr[0]:.3%}, PC2 {evr[1]:.3%}",
           f"groups (GMM on PC1): left {np.sum(g == 0)}, right {np.sum(g == 1)}", ""]

    # --- composition
    for col in ("State", "Sex", "Site"):
        t = pd.crosstab(info[col], info.Group)
        p = chi2_contingency(t)[1]
        t.to_csv(OUT / f"group_vs_{col.lower()}.csv")
        if col != "Site":
            rep += [f"Group x {col}: chi2 p = {p:.3g}", t.to_string(), ""]
        else:
            rep += [f"Group x Site ({t.shape[0]} sites): chi2 p = {p:.3g}", ""]

    # --- batch-like per-subject QC on ALL SNPs of the .bed (chunked)
    miss = np.zeros(len(rows)); het = np.zeros(len(rows)); nobs = np.zeros(len(rows))
    nb = len(pd.read_csv(f"{PFX}.bim", sep=r"\s+", header=None, usecols=[0]))
    for i in range(0, nb, 20000):
        Gq = U.read_bed(PFX, range(i, min(i + 20000, nb)), n_all)[rows]
        m = np.isnan(Gq)
        miss += m.sum(1); nobs += (~m).sum(1); het += (Gq == 1).sum(1)
    info["Missing_rate"] = miss / nb; info["Het_rate"] = het / nobs
    for c in ("Missing_rate", "Het_rate"):
        a, b = info.loc[info.Group == "left", c], info.loc[info.Group == "right", c]
        rep.append(f"{c}: left median {a.median():.5f}, right median {b.median():.5f}, "
                   f"Mann-Whitney p = {mannwhitneyu(a, b).pvalue:.3g}")
    rep.append("")
    info.to_csv(OUT / "subjects.csv", index=False)

    # --- PC1 drivers
    snp2reg = {}
    for r, ii in r2s.items():
        for i in ii:
            snp2reg.setdefault(i, []).append(r)
    L = pd.DataFrame({"SNP": bim.loc[snp, "snp"].values if "snp" in bim else bim.iloc[snp, 1].values,
                      "chr": bim.loc[snp, "chr"].values, "bp": bim.loc[snp, "bp"].values,
                      "Regions": [";".join(snp2reg.get(i, [])) for i in snp],
                      "PC1_loading2": Vt[0] ** 2})
    fl, fr = X[g == 0].mean(0), X[g == 1].mean(0)
    L["Carrier_freq_left"], L["Carrier_freq_right"], L["Diff"] = fl, fr, fr - fl
    L["Fisher_p"] = [fisher_exact([[X[g == 0, j].sum(), (g == 0).sum() - X[g == 0, j].sum()],
                                   [X[g == 1, j].sum(), (g == 1).sum() - X[g == 1, j].sum()]])[1]
                     for j in range(X.shape[1])]
    L = L.sort_values("PC1_loading2", ascending=False)
    L.to_csv(OUT / "snps_pc1.csv", index=False)
    bychr = L.groupby("chr").agg(SNPs=("SNP", "size"), PC1_share=("PC1_loading2", "sum")).sort_values("PC1_share", ascending=False)
    bychr.to_csv(OUT / "pc1_by_chr.csv")
    rep += ["PC1 share by chromosome (sum of squared loadings, top 8):", bychr.head(8).round(4).to_string(), ""]
    top = L.head(50)
    win = top.groupby("chr").bp.agg(["min", "max", "size"])
    rep += ["top-50 PC1 SNPs by chromosome (span):", win.to_string(), ""]
    rr = L.assign(Region=L.Regions.str.split(";")).explode("Region").groupby("Region") \
        .agg(chr=("chr", "first"), start=("bp", "min"), end=("bp", "max"), SNPs=("SNP", "size"),
             PC1_share=("PC1_loading2", "sum"), max_absdiff=("Diff", lambda d: d.abs().max())) \
        .sort_values("PC1_share", ascending=False)
    rr.to_csv(OUT / "regions_pc1.csv")
    rep += ["top regions by PC1 share:", rr.head(15).round(4).to_string(), "",
            f"SNPs with Fisher p < 1e-10 between groups: {(L.Fisher_p < 1e-10).sum()}; "
            f"|carrier-freq diff| > 0.5: {(L.Diff.abs() > 0.5).sum()}", ""]

    # --- how few SNPs discriminate the groups
    cv = StratifiedKFold(5, shuffle=True, random_state=1)
    order = np.argsort(-Vt[0] ** 2)
    rep.append("cross-validated accuracy predicting left/right (logistic regression; SNPs ranked by PC1 loading on all "
               "subjects, so mildly optimistic):")
    for k in (1, 3, 10, 30, 100):
        acc = cross_val_score(LogisticRegression(max_iter=2000), X[:, order[:k]], g, cv=cv).mean()
        rep.append(f"  top {k:>3} SNPs: {acc:.3f}")
    reg0 = rr.index[0]
    cols = [j for j in range(X.shape[1]) if reg0 in L.sort_index().Regions.values[j].split(";")]
    acc = cross_val_score(LogisticRegression(max_iter=2000), X[:, cols], g, cv=cv).mean()
    rep += [f"  only region {reg0} ({len(cols)} SNPs): {acc:.3f}",
            f"  majority class baseline: {max(g.mean(), 1 - g.mean()):.3f}", ""]

    # --- same split in other feature sets?
    rng = np.random.default_rng(1); pool = sorted(r2s)
    others = {"SKAT": sorted(U.selected_regions(U.SELECTIONS[("SKAT", REGION)]) & set(r2s)),
              "Random@DiGAS d0": list(rng.choice(pool, len(genes), replace=False)), "AllRegions": pool}
    top_regions = set(rr.index[:10])
    rep.append("same subjects split along PC1 of other sets? (ARI between DiGAS groups and a GMM split of that set's PC1)")
    for name, gs in others.items():
        ii = sorted({i for r in gs for i in r2s[r]})
        Zo = U.pca_scores(PFX, ii, n_all, rows, one, n_comp=5)["all"][0]
        go = split_pc1(Zo)
        rep.append(f"  {name}: ARI {adjusted_rand_score(g, go):.3f}; |corr(PC1, DiGAS PC1)| "
                   f"{abs(np.corrcoef(Zo[:, 0], Z[:, 0])[0, 1]):.3f}; contains {len(top_regions & set(gs))}/10 "
                   f"of DiGAS's top PC1 regions")
    open(OUT / "report.txt", "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
