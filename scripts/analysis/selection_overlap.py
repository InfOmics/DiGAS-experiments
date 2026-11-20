#!/usr/bin/env python3
"""Which regions do DiGAS and SKAT select in common? Quantitative part (whole-cohort selections).

Cases: exon_5k in ADNI1/ADNI2/ADNI3, and ADNI1 exon / exon_protein_coding.
Per case and contrast (CN-AD, CN-MCI, and their union = the set used downstream):
  1. overlap of the selected sets (p <= 0.05): counts, Jaccard, overlap coefficient, fold over chance,
     hypergeometric p, and the expected overlap when both tools' selection rates depend on region size
     (stratified by SNP count, deciles) -> does agreement exceed what a shared size bias produces?
  2. selection rate by region size (SNP-count bins) for each tool
  3. agreement over ALL tested regions: Spearman of the two tools' p-values, top-k overlap vs chance
  4. cross-cohort: each tool's exon_5k selection in ADNI1 vs ADNI2 vs ADNI3, and the common (DiGAS & SKAT) set
  5. link to the unsupervised analysis (CN+AD subjects): 10-NN excess in 50 PCs and k-means ARI for the
     Common, DiGAS-only and SKAT-only region sets, each with a size-matched random reference (5 draws)
No label permutations are used: references are analytic (hypergeometric) or label-free random region sets.

  .venv/bin/python scripts/analysis/selection_overlap.py
Outputs: results/selection_overlap/
"""
import importlib.util, os
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("U", ROOT / "scripts/analysis/unsupervised_tools.py")
U = importlib.util.module_from_spec(spec); spec.loader.exec_module(U)
OUT = ROOT / "results/selection_overlap"; OUT.mkdir(parents=True, exist_ok=True)
SEL = ROOT / "results/tool_comparison/selection_all"
CASES = {  # (dataset, region) -> (DiGAS prefix, SKAT prefix)
    ("ADNI1", "exon_5k"): (ROOT / "results/fold_pair_scatter_leaky/digas_all_subjects/digas_all",
                           ROOT / "results/skat_leaky/skat_all_exon_5k"),
    ("ADNI2", "exon_5k"): (SEL / "DiGAS_exon_5k_ADNI2", SEL / "SKAT_exon_5k_ADNI2"),
    ("ADNI3", "exon_5k"): (SEL / "DiGAS_exon_5k_ADNI3", SEL / "SKAT_exon_5k_ADNI3"),
    ("ADNI1", "exon"): (SEL / "DiGAS_exon", SEL / "SKAT_exon"),
    ("ADNI1", "exon_protein_coding"): (SEL / "DiGAS_exon_protein_coding", SEL / "SKAT_exon_protein_coding"),
}
CONTRASTS = ["CN-AD", "CN-MCI"]


def load_all(prefix):
    """all tested regions: Gene x contrast -> p (contrast named CN-AD / CN-MCI / AD-MCI regardless of order)."""
    d = pd.read_csv(f"{prefix}_all.txt", sep=r"\s+", dtype=str)
    d["P"] = pd.to_numeric(d.Pvalue, errors="coerce")
    d["Contrast"] = [("CN-" + (b if a == "CN" else a)) if "CN" in (a, b) else "AD-MCI" for a, b in zip(d.Cat1, d.Cat2)]
    p = d.pivot_table(index="Gene", columns="Contrast", values="P", aggfunc="first")
    if "N_markers_tested" in d:
        p["N_snps"] = pd.to_numeric(d.drop_duplicates("Gene").set_index("Gene").N_markers_tested, errors="coerce")
    return p


def overlap_row(A, B, N, size=None):
    k, a, b = len(A & B), len(A), len(B)
    exp = a * b / N
    r = dict(N_tested=N, DiGAS=a, SKAT=b, Common=k, Expected=exp, Fold=k / exp if exp else np.nan,
             Hypergeom_p=hypergeom.sf(k - 1, N, a, b), Jaccard=k / len(A | B) if A | B else np.nan,
             Overlap_coef=k / min(a, b) if min(a, b) else np.nan)
    if size is not None:   # expected overlap if selection depends only on region size (independent within bins)
        bins = pd.qcut(size.rank(method="first"), 10, labels=False)
        e = 0.0
        for _, g in size.groupby(bins):
            idx = set(g.index); n = len(idx)
            e += len(A & idx) * len(B & idx) / n
        r["Expected_size_strat"] = e; r["Fold_size_strat"] = k / e if e else np.nan
    return r


def main():
    rows, rates, conc, topk, sets = [], [], [], [], {}
    for (ds, reg), (pd_, ps_) in CASES.items():
        D, S = load_all(pd_), load_all(ps_)
        common_idx = D.index.intersection(S.index)
        size = S.loc[common_idx, "N_snps"]
        union = {"DiGAS": set(), "SKAT": set()}
        for c in CONTRASTS:
            ok = D.loc[common_idx, c].notna() & S.loc[common_idx, c].notna()
            u = common_idx[ok.values]
            A = set(u[D.loc[u, c].values <= 0.05]); B = set(u[S.loc[u, c].values <= 0.05])
            union["DiGAS"] |= A; union["SKAT"] |= B
            rows.append(dict(Dataset=ds, Region=reg, Contrast=c, **overlap_row(A, B, len(u), size.loc[u])))
            pdg, psk = D.loc[u, c].values, S.loc[u, c].values
            rho = spearmanr(pdg, psk)
            conc.append(dict(Dataset=ds, Region=reg, Contrast=c, N=len(u), Spearman_p=rho.correlation,
                             Spearman_pval=rho.pvalue))
            rd, rs = pd.Series(pdg, u).rank(method="first"), pd.Series(psk, u).rank(method="first")
            for k in (50, 100, 250, 500, 1000, 2000):
                if k < len(u):
                    kk = len(set(rd[rd <= k].index) & set(rs[rs <= k].index))
                    topk.append(dict(Dataset=ds, Region=reg, Contrast=c, Top_k=k, Common=kk, Expected=k * k / len(u),
                                     Fold=kk / (k * k / len(u))))
            for t, sel in (("DiGAS", A), ("SKAT", B)):
                sb = pd.cut(size.loc[u], [0, 1, 2, 5, 10, 20, 50, 1e9], labels=["1", "2", "3-5", "6-10", "11-20", "21-50", ">50"])
                fr = pd.Series(u.isin(list(sel)), u).groupby(sb.values).agg(["mean", "size"])
                for b_, r_ in fr.iterrows():
                    rates.append(dict(Dataset=ds, Region=reg, Contrast=c, Tool=t, SNP_bin=b_, Rate=r_["mean"], N=r_["size"]))
        A, B = union["DiGAS"], union["SKAT"]
        rows.append(dict(Dataset=ds, Region=reg, Contrast="union", **overlap_row(A, B, len(common_idx), size)))
        sets[(ds, reg)] = dict(DiGAS=A, SKAT=B, Common=A & B, universe=set(common_idx))
        pd.DataFrame({"Gene": sorted(A | B), "DiGAS": [g in A for g in sorted(A | B)],
                      "SKAT": [g in B for g in sorted(A | B)]}).to_csv(OUT / f"selected_{ds}_{reg}.csv", index=False)
        # per-region table for the common set (for the qualitative part)
        cm = sorted(A & B)
        t = pd.DataFrame(index=cm)
        for c in CONTRASTS:
            t[f"DiGAS_p_{c}"] = D.loc[cm, c]; t[f"SKAT_p_{c}"] = S.loc[cm, c]
        t["N_snps"] = size.loc[cm]
        t.index.name = "Gene"; t.to_csv(OUT / f"common_{ds}_{reg}.csv")

    ov = pd.DataFrame(rows); ov.to_csv(OUT / "overlap.csv", index=False)
    pd.DataFrame(rates).to_csv(OUT / "selection_rate_by_size.csv", index=False)
    pd.DataFrame(conc).to_csv(OUT / "pvalue_concordance.csv", index=False)
    pd.DataFrame(topk).to_csv(OUT / "topk_overlap.csv", index=False)

    # cross-cohort (exon_5k)
    xc = []
    cohorts = ["ADNI1", "ADNI2", "ADNI3"]
    for i, a in enumerate(cohorts):
        for b in cohorts[i + 1:]:
            sa, sb = sets[(a, "exon_5k")], sets[(b, "exon_5k")]
            uni = sa["universe"] & sb["universe"]
            for t in ("DiGAS", "SKAT", "Common"):
                xc.append(dict(Pair=f"{a}-{b}", Set=t, **overlap_row(sa[t] & uni, sb[t] & uni, len(uni))))
    xc = pd.DataFrame(xc).rename(columns={"DiGAS": "N_first", "SKAT": "N_second"}); xc.to_csv(OUT / "cross_cohort.csv", index=False)
    s3 = [sets[(c, "exon_5k")] for c in cohorts]
    uni = set.intersection(*[s["universe"] for s in s3])
    rep3 = {t: sorted(set.intersection(*[s[t] for s in s3]) & uni) for t in ("DiGAS", "SKAT", "Common")}
    pd.Series({t: ";".join(v) for t, v in rep3.items()}).to_csv(OUT / "in_all_three_cohorts.csv", header=["Genes"])

    # unsupervised link
    un = []
    for ds in cohorts:
        s = sets[(ds, "exon_5k")]
        un += unsup(ds, {"Common": s["Common"], "DiGAS-only": s["DiGAS"] - s["SKAT"], "SKAT-only": s["SKAT"] - s["DiGAS"],
                         "DiGAS (all)": s["DiGAS"], "SKAT (all)": s["SKAT"]})
    un = pd.DataFrame(un); un.to_csv(OUT / "unsupervised_subsets.csv", index=False)

    pd.set_option("display.width", 220)
    rep = ["== overlap of selected regions", ov.round(4).to_string(index=False), "",
           "== Spearman of p-values over all tested regions", pd.DataFrame(conc).round(4).to_string(index=False), "",
           "== top-k overlap (CN-AD)", pd.DataFrame(topk).query("Contrast=='CN-AD'").round(2).to_string(index=False), "",
           "== cross-cohort (exon_5k)", xc.round(4).to_string(index=False), "",
           "== regions selected in all three cohorts: " + ", ".join(f"{t} {len(v)}" for t, v in rep3.items()), "",
           "== unsupervised on subsets (CN+AD)", un.round(3).to_string(index=False)]
    open(OUT / "report.txt", "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


def unsup(ds, groups, draws=5, seed=1):
    reg = "exon_5k"; pfx = ROOT / f"{ds}_QC"
    subj, bim, regions = U.C.load_data(ROOT, f"{ds}_QC", ROOT / U.C.DATASETS[f"{ds}_QC"], ROOT / "region" / f"{reg}.txt")
    fam = pd.read_csv(f"{pfx}.fam", sep=r"\s+", header=None, usecols=[1], names=["IID"], dtype=str)
    lab = subj.set_index(subj.IID.astype(str)).OriginalGroup
    ra = np.where(fam.IID.isin(lab.index))[0]; st = lab.loc[fam.IID.iloc[ra]].values
    rows = ra[np.isin(st, ["CN", "AD"])]; y = lab.loc[fam.IID.iloc[rows]].values
    r2s = U.region_snps(regions, bim.reset_index(drop=True)); pool = sorted(r2s)
    rng = np.random.default_rng(seed); one = {"x": np.ones(len(rows), bool)}
    out = []

    def score(genes):
        idx = sorted({i for g in genes for i in r2s[g]})
        Z = U.pca_scores(pfx, idx, len(fam), rows, one)["x"][0]
        o, e = U.knn_agreement(Z, y)
        ari = max(adjusted_rand_score(y, KMeans(k, n_init=20, random_state=seed).fit_predict(Z)) for k in (2, 3))
        return o - e, ari, len(idx)

    for name, g in groups.items():
        g = sorted(set(g) & set(r2s))
        if len(g) < 5:
            out.append(dict(Dataset=ds, Set=name, Regions=len(g))); continue
        k, a, n = score(g)
        rk = [score(rng.choice(pool, len(g), replace=False))[:2] for _ in range(draws)]
        out.append(dict(Dataset=ds, Set=name, Regions=len(g), SNPs=n, kNN10_excess=k, Best_ARI=a,
                        Random_kNN_mean=np.mean([r[0] for r in rk]), Random_kNN_max=np.max([r[0] for r in rk]),
                        Random_ARI_max=np.max([r[1] for r in rk])))
    return out


if __name__ == "__main__":
    main()
