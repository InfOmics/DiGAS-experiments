"""Unsupervised, model-free variability of ADNI1_QC by diagnosis group (CN / MCI / AD).

No classifier and no shuffled labels: every reference (null) distribution is analytic sampling theory.

Stage 0  structure (PLINK, autosomes, MAF>=0.05, LD-pruned 50/5/0.2):
         PCA (10 PCs) vs group, pairwise kinship (PI_HAT) within vs between groups, heterozygosity F by group, sex by group.
Stage 1  per SNP x group (autosomes): allele frequency (cohort minor allele), carrier frequency, observed/expected
         heterozygosity, missingness, HWE exact p (PLINK --hardy per group), binomial SE of the allele frequency and
         a 200-replicate within-group subject bootstrap SD.
Stage 2  between groups, per SNP: allelic chi2 for CN/MCI/AD (df 2) and the pairs CN-AD, CN-MCI, MCI-AD, CN-(AD+MCI)
         (df 1); Hudson Fst; differential missingness. Whole cohort: lambda_GC, Storey pi1 (share of SNPs whose
         frequency truly differs), observed vs expected tail counts, genome-wide Fst with chromosome jackknife SE.
         Done on all autosomal SNPs and on the LD-pruned set (approximately independent SNPs).
Stage 3  what the dataset can resolve: per SNP, the smallest allelic odds ratio detectable with 80% power
         (alpha 0.05 and 5e-8) for each group pair, given its frequency and the group sizes; sample size needed
         for typical AD-locus effects.
Stage 4  gene level (exon_5k, DiGAS's unit): per-subject Hill value (cohort minor allele, as in
         region_feature_pipeline.py), group Hill, bootstrap SE of the CN-AD and CN-MCI log1p difference (|.| = DiGAS FC),
         z = difference / SE, and how SE and the top-FC ranks depend on the number of SNPs per gene.

Outputs: results/snp_variability/
Usage: python scripts/analysis/snp_variability.py
"""
import os, sys, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

ROOT = Path(HERE).parents[1]
OUT = ROOT / "results" / "snp_variability"; TMP = OUT / "plink"
GROUPS = ("CN", "MCI", "AD")
COL = {"CN": "#2a78d6", "MCI": "#eb6834", "AD": "#1baf7a"}          # validated categorical slots 1-3
INK, MUTED = "#0b0b0b", "#52514e"
PAIRS = [("CN", "AD"), ("CN", "MCI"), ("MCI", "AD"), ("CN", "AD+MCI")]
NBOOT, CHUNK, SEED = 200, 40000, 0
KNOWN = {"ENSG00000130203": "APOE", "ENSG00000130204": "TOMM40", "ENSG00000130208": "APOC1",
         "ENSG00000130202": "NECTIN2", "ENSG00000267282": "lncRNA (chr19, near TOMM40)"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#ecebe7", "grid.linewidth": 0.6, "figure.dpi": 150})


def plink(*args):
    subprocess.run([str(ROOT / "plink"), "--bfile", str(ROOT / "ADNI1_QC"), "--allow-no-sex", *args],
                   check=True, stdout=subprocess.DEVNULL)


def savefig(fig, name):
    fig.tight_layout(); fig.savefig(OUT / f"{name}.png", dpi=160); fig.savefig(OUT / f"{name}.pdf"); plt.close(fig)


def load_meta():
    fam = pd.read_csv(ROOT / "ADNI1_QC.fam", sep=r"\s+", header=None, names=["fid", "iid", "pat", "mat", "sex", "ph"])
    ph = pd.read_csv(ROOT / "pheno/adni1pheno.txt", sep="\t", header=None, names=["iid", "group"])
    fam = fam.merge(ph, on="iid", how="left"); assert fam.group.notna().all()
    bim = pd.read_csv(ROOT / "ADNI1_QC.bim", sep="\t", header=None, names=["chr", "snp", "cm", "bp", "a1", "a2"])
    return fam, bim


# ---------------- stage 0: structure ----------------
def structure(fam):
    TMP.mkdir(parents=True, exist_ok=True)
    for g in GROUPS:
        fam.loc[fam.group == g, ["fid", "iid"]].to_csv(TMP / f"keep_{g}.txt", sep=" ", header=False, index=False)
        if not (TMP / f"hardy_{g}.hwe").exists():
            plink("--autosome", "--keep", str(TMP / f"keep_{g}.txt"), "--hardy", "--out", str(TMP / f"hardy_{g}"))
    if not (TMP / "prune.prune.in").exists():
        plink("--autosome", "--maf", "0.05", "--indep-pairwise", "50", "5", "0.2", "--out", str(TMP / "prune"))
    pr = ["--autosome", "--extract", str(TMP / "prune.prune.in")]
    if not (TMP / "pca.eigenvec").exists(): plink(*pr, "--pca", "10", "--out", str(TMP / "pca"))
    if not (TMP / "genome.genome").exists(): plink(*pr, "--genome", "--out", str(TMP / "genome"))
    if not (TMP / "het.het").exists(): plink(*pr, "--het", "--out", str(TMP / "het"))

    rows = []
    pc = pd.read_csv(TMP / "pca.eigenvec", sep=r"\s+", header=None,
                     names=["fid", "iid"] + [f"PC{i}" for i in range(1, 11)]).merge(fam[["iid", "group"]], on="iid")
    ev = np.loadtxt(TMP / "pca.eigenval")
    for i in range(1, 11):
        x = [pc.loc[pc.group == g, f"PC{i}"].values for g in GROUPS]
        rows.append(dict(test=f"PC{i} by group (Kruskal-Wallis)", eigenvalue=ev[i - 1],
                         stat=stats.kruskal(*x).statistic, p=stats.kruskal(*x).pvalue,
                         **{f"median_{g}": np.median(v) for g, v in zip(GROUPS, x)}))
    het = pd.read_csv(TMP / "het.het", sep=r"\s+").rename(columns={"IID": "iid"}).merge(fam[["iid", "group"]], on="iid")
    x = [het.loc[het.group == g, "F"].values for g in GROUPS]
    rows.append(dict(test="heterozygosity F by group (Kruskal-Wallis)", stat=stats.kruskal(*x).statistic,
                     p=stats.kruskal(*x).pvalue, **{f"median_{g}": np.median(v) for g, v in zip(GROUPS, x)}))
    ct = pd.crosstab(fam.group, fam.sex).loc[list(GROUPS)]
    c2 = stats.chi2_contingency(ct)
    rows.append(dict(test="sex by group (chi2)", stat=c2[0], p=c2[1],
                     **{f"median_{g}": ct.loc[g, 2] / ct.loc[g].sum() for g in GROUPS}))   # share female (code 2)
    gm = pd.read_csv(TMP / "genome.genome", sep=r"\s+", usecols=["IID1", "IID2", "PI_HAT"])
    grp = fam.set_index("iid").group
    gm["g1"] = gm.IID1.map(grp).values; gm["g2"] = gm.IID2.map(grp).values
    gm["kind"] = np.where(gm.g1 == gm.g2, "within " + gm.g1, "between")
    kin = gm.groupby("kind").PI_HAT.agg(pairs="size", mean="mean", median="median",
                                        related_0_1875=lambda v: int((v > 0.1875).sum()),
                                        related_0_1=lambda v: int((v > 0.1).sum())).reset_index()
    kin.to_csv(OUT / "kinship_within_between.csv", index=False)
    gm[gm.PI_HAT > 0.1].to_csv(TMP / "related_pairs.csv", index=False)   # subject IDs: not in results/
    S = pd.DataFrame(rows); S.to_csv(OUT / "structure_tests.csv", index=False)

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    for k, (a, b) in enumerate([(1, 2), (3, 4)]):
        for g in GROUPS:
            d = pc[pc.group == g]
            ax[k].scatter(d[f"PC{a}"], d[f"PC{b}"], s=9, c=COL[g], label=f"{g} (n={len(d)})", alpha=.75,
                          edgecolors="white", linewidths=.4)
        ax[k].set_xlabel(f"PC{a}"); ax[k].set_ylabel(f"PC{b}")
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title("Ancestry PCs by group (LD-pruned autosomes)", loc="left", fontsize=9, color=INK)
    parts = ax[2].violinplot(x, showmedians=True, widths=.7)
    for b_, g in zip(parts["bodies"], GROUPS): b_.set_facecolor(COL[g]); b_.set_alpha(.6)
    ax[2].set_xticks([1, 2, 3], GROUPS); ax[2].set_ylabel("heterozygosity F (inbreeding)")
    ax[2].set_title("Per-subject heterozygosity", loc="left", fontsize=9, color=INK)
    savefig(fig, "structure_pca_het")
    return S, kin


# ---------------- stage 1-3: per SNP ----------------
def decode(raw, nb, n, lo, hi):
    M = raw[lo:hi]
    c = np.stack([(M >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(hi - lo, nb * 4)[:, :n]
    dos = np.select([c == 0, c == 2, c == 3], [2, 1, 0], -1).astype(np.int8)   # copies of A1, -1 = missing
    return dos


def chi2_allelic(c, n2):
    """c: (S, k) counts of allele, n2: (S, k) total alleles -> chi2 (df k-1), vectorised Pearson on the 2 x k table."""
    tot = n2.sum(1, keepdims=True); p = c.sum(1, keepdims=True) / tot
    e1, e0 = n2 * p, n2 * (1 - p)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = ((c - e1) ** 2 / e1 + ((n2 - c) - e0) ** 2 / e0).sum(1)
    return np.where(np.isfinite(x), x, np.nan)


def hudson(p1, p2, n1, n2):
    """Hudson Fst numerator/denominator (Bhatia et al. 2013); n = allele counts."""
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def min_detectable_or(p0, n0, n1, alpha, power=0.8):
    """Smallest OR>1 for the allelic test (n0, n1 = subjects) at reference allele freq p0; bisection on log OR."""
    za, zb = stats.norm.isf(alpha / 2), stats.norm.ppf(power)
    lo, hi = np.zeros_like(p0), np.full_like(p0, np.log(50.0))
    for _ in range(50):
        m = (lo + hi) / 2; o = np.exp(m); p1 = o * p0 / (1 - p0 + o * p0)
        pb = (n0 * p0 + n1 * p1) / (n0 + n1)
        se0 = np.sqrt(pb * (1 - pb) * (1 / (2 * n0) + 1 / (2 * n1)))
        se1 = np.sqrt(p0 * (1 - p0) / (2 * n0) + p1 * (1 - p1) / (2 * n1))
        ok = np.abs(p1 - p0) >= za * se0 + zb * se1
        hi = np.where(ok, m, hi); lo = np.where(ok, lo, m)
    return np.exp(hi)


def n_needed(p0, OR, alpha, power=0.8):
    """Subjects per group (balanced) for the allelic test."""
    p1 = OR * p0 / (1 - p0 + OR * p0); pb = (p0 + p1) / 2
    za, zb = stats.norm.isf(alpha / 2), stats.norm.ppf(power)
    a = za * np.sqrt(2 * pb * (1 - pb)) + zb * np.sqrt(p0 * (1 - p0) + p1 * (1 - p1))
    return int(np.ceil(a ** 2 / (p1 - p0) ** 2 / 2))   # alleles -> subjects


def storey_pi1(p, lam=0.5):
    p = p[np.isfinite(p)]
    return 1 - min(1.0, (p > lam).mean() / (1 - lam))


def per_snp(fam, bim):
    auto = np.where(bim.chr.between(1, 22))[0]
    n = len(fam); nb = (n + 3) // 4
    raw = np.fromfile(ROOT / "ADNI1_QC.bed", dtype=np.uint8); assert raw[:3].tolist() == [0x6c, 0x1b, 1]
    raw = raw[3:].reshape(len(bim), nb)
    idx = {g: np.where(fam.group.values == g)[0] for g in GROUPS}
    rng = np.random.default_rng(SEED)
    W = {g: rng.multinomial(len(i), np.full(len(i), 1 / len(i)), size=NBOOT).astype(np.float32) for g, i in idx.items()}
    out = []
    t0 = time.time()
    for lo in range(0, len(auto), CHUNK):
        sel = auto[lo:lo + CHUNK]; dos = decode(raw, nb, n, sel[0], sel[-1] + 1)[sel - sel[0]]
        r = {"snp": bim.snp.values[sel], "chr": bim.chr.values[sel], "bp": bim.bp.values[sel]}
        a_all = np.where(dos >= 0, dos, 0).sum(1); n_all = (dos >= 0).sum(1)
        flip = a_all / np.maximum(2 * n_all, 1) > 0.5          # count the cohort minor allele
        for g, i in idx.items():
            d = dos[:, i]; obs = d >= 0
            a = np.where(obs, d, 0).astype(np.float32); a = np.where(flip[:, None] & obs, 2 - a, a)
            no = obs.sum(1); c = a.sum(1); p = c / (2 * no)
            r[f"n_{g}"] = no; r[f"ac_{g}"] = c; r[f"maf_{g}"] = p
            r[f"carrier_{g}"] = ((a >= 1) & obs).sum(1) / no
            r[f"het_obs_{g}"] = ((a == 1) & obs).sum(1) / no; r[f"het_exp_{g}"] = 2 * p * (1 - p)
            r[f"miss_{g}"] = 1 - no / len(i)
            r[f"se_binom_{g}"] = np.sqrt(p * (1 - p) / (2 * no))
            pb = (W[g] @ a.T) / (2 * (W[g] @ obs.T.astype(np.float32)))           # NBOOT x S
            r[f"se_boot_{g}"] = pb.std(0, ddof=1)
        out.append(pd.DataFrame(r))
        print(f"  SNPs {lo + len(sel):>7}/{len(auto)}  {time.time() - t0:.0f}s", flush=True)
    D = pd.concat(out, ignore_index=True)
    for g in GROUPS:                                                               # HWE per group from PLINK
        h = pd.read_csv(TMP / f"hardy_{g}.hwe", sep=r"\s+", usecols=["SNP", "P"]).rename(columns={"P": f"hwe_p_{g}"})
        D = D.merge(h, left_on="snp", right_on="SNP", how="left").drop(columns="SNP")
    D["ac_AD+MCI"] = D.ac_AD + D.ac_MCI; D["n_AD+MCI"] = D.n_AD + D.n_MCI
    D["maf_AD+MCI"] = D["ac_AD+MCI"] / (2 * D["n_AD+MCI"])
    C = np.c_[[D[f"ac_{g}"] for g in GROUPS]].T; N2 = 2 * np.c_[[D[f"n_{g}"] for g in GROUPS]].T
    D["chi2_3groups"] = chi2_allelic(C, N2); D["p_3groups"] = stats.chi2.sf(D.chi2_3groups, 2)
    for a, b in PAIRS:
        k = f"{a}_vs_{b}"
        x = chi2_allelic(np.c_[D[f"ac_{a}"], D[f"ac_{b}"]], 2 * np.c_[D[f"n_{a}"], D[f"n_{b}"]])
        D[f"chi2_{k}"] = x; D[f"p_{k}"] = stats.chi2.sf(x, 1)
        D[f"diff_{k}"] = D[f"maf_{b}"] - D[f"maf_{a}"]
        num, den = hudson(D[f"maf_{a}"], D[f"maf_{b}"], 2 * D[f"n_{a}"], 2 * D[f"n_{b}"])
        D[f"fst_num_{k}"] = num; D[f"fst_den_{k}"] = den
    M = np.c_[[D[f"miss_{g}"] * (fam.group == g).sum() for g in GROUPS]].T          # missing counts
    D["p_diffmiss"] = stats.chi2.sf(chi2_allelic(M, np.c_[[np.full(len(D), (fam.group == g).sum()) for g in GROUPS]].T), 2)
    return D


def cohort_summary(D, pruned):
    rows = []
    for lab, X in (("all autosomal SNPs", D), ("LD-pruned (MAF>=0.05, r2<0.2)", D[D.snp.isin(pruned)])):
        for a, b in PAIRS + [("3", "groups")]:
            k = "3groups" if a == "3" else f"{a}_vs_{b}"; p = X[f"p_{k}"].values; p = p[np.isfinite(p)]
            x = stats.chi2.isf(p, 2 if k == "3groups" else 1)
            r = dict(snp_set=lab, comparison=k.replace("_vs_", " vs "), snps=len(p),
                     lambda_gc=np.median(x) / stats.chi2.ppf(.5, 2 if k == "3groups" else 1), pi1_storey=storey_pi1(p))
            for t in (1e-3, 1e-4, 1e-5, 1e-6, 5e-8):
                r[f"obs_p<{t:g}"] = int((p < t).sum()); r[f"exp_p<{t:g}"] = round(t * len(p), 2)
            if k != "3groups":                                   # genome-wide Hudson Fst, chromosome jackknife
                num, den, ch = X[f"fst_num_{k}"].values, X[f"fst_den_{k}"].values, X.chr.values
                ok = np.isfinite(num) & np.isfinite(den) & (den > 0); num, den, ch = num[ok], den[ok], ch[ok]
                f = num.sum() / den.sum()
                jk = np.array([(num[ch != c].sum() / den[ch != c].sum()) for c in range(1, 23)])
                se = np.sqrt((len(jk) - 1) / len(jk) * ((jk - jk.mean()) ** 2).sum())
                r.update(fst_hudson=f, fst_se_jackknife=se, fst_z=f / se)
            rows.append(r)
    for g in GROUPS:                                             # bootstrap SD vs binomial SE, HWE
        ok = D[f"maf_{g}"] > 0.01
        rows.append(dict(snp_set="all autosomal SNPs", comparison=f"within {g}",
                         snps=int(ok.sum()), boot_over_binom_SE_median=np.median(D.loc[ok, f"se_boot_{g}"] / D.loc[ok, f"se_binom_{g}"]),
                         se_binom_median=D.loc[ok, f"se_binom_{g}"].median(),
                         **{f"obs_hwe_p<{t:g}": int((D[f"hwe_p_{g}"] < t).sum()) for t in (1e-3, 1e-6)},
                         **{f"exp_hwe_p<{t:g}": round(t * D[f"hwe_p_{g}"].notna().sum(), 1) for t in (1e-3, 1e-6)}))
    p = D.p_diffmiss.values; p = p[np.isfinite(p)]
    rows.append(dict(snp_set="all autosomal SNPs", comparison="differential missingness CN/MCI/AD", snps=len(p),
                     pi1_storey=storey_pi1(p), **{f"obs_p<{t:g}": int((p < t).sum()) for t in (1e-3, 1e-6)},
                     **{f"exp_p<{t:g}": round(t * len(p), 2) for t in (1e-3, 1e-6)}))
    return pd.DataFrame(rows)


def power_tables(D, fam):
    ng = fam.group.value_counts(); ng["AD+MCI"] = ng.AD + ng.MCI
    rows, need = [], []
    bins = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
    for a, b in PAIRS:
        p0 = D[f"maf_{a}"].clip(1e-4, .5).values
        for alpha in (0.05, 5e-8):
            mo = min_detectable_or(p0, ng[a], ng[b], alpha)
            D[f"minOR_{a}_vs_{b}_{alpha:g}"] = mo
            cut = pd.cut(p0, bins)
            for iv, v in pd.Series(mo).groupby(cut, observed=True):
                rows.append(dict(comparison=f"{a} vs {b}", n=f"{ng[a]} vs {ng[b]}", alpha=alpha, maf_bin=str(iv),
                                 snps=len(v), min_detectable_OR_median=v.median()))
    for OR in (1.1, 1.15, 1.2, 1.5, 3.0):
        for p0 in (0.05, 0.1, 0.2, 0.3, 0.4):
            need.append(dict(OR=OR, maf=p0, n_per_group_alpha_0_05=n_needed(p0, OR, .05),
                             n_per_group_alpha_5e_8=n_needed(p0, OR, 5e-8)))
    return pd.DataFrame(rows), pd.DataFrame(need)


def gene_map(top, region="gene"):
    reg = pd.read_csv(ROOT / f"region/{region}.txt", sep="\t", header=None, names=["g", "s", "e"])
    reg["chr"] = reg.s.str.split(":").str[0].astype(int)
    reg["s"] = reg.s.str.split(":").str[1].astype(int); reg["e"] = reg.e.str.split(":").str[1].astype(int)
    lab = []
    for c, bp in zip(top.chr, top.bp):
        r = reg[(reg.chr == c) & (reg.s - 20000 <= bp) & (reg.e + 20000 >= bp)]
        lab.append(";".join(KNOWN.get(g.split(".")[0], g.split(".")[0]) for g in r.g.unique()[:4]))
    return lab



# ---------------- per-subject region values (same code as region_feature_pipeline.py of the working repository) ----------------
def carriers(D, tr_mask):
    """SNP x subject carrier matrix (float32), counted allele = minor allele among the training subjects."""
    hp = D["hap"][:, tr_mask]
    h1 = D["hom1"][:, tr_mask]; h2 = D["hom2"][:, tr_mask]; he = D["het"][:, tr_mask]
    hd = (he & ~hp).sum(1)
    a1 = 2 * h1.sum(1) + hd - (h1 & hp).sum(1); a2 = 2 * h2.sum(1) + hd - (h2 & hp).sum(1)
    swap = (a1 > a2)[:, None]
    return np.where(swap, D["hom2"] | D["het"], D["hom1"] | D["het"]).astype(np.float32)


def region_values(D, tr_mask):
    """regions x subjects: per-subject Hill score."""
    return np.asarray(D["W"].T @ carriers(D, tr_mask), dtype=np.float32)


# ---------------- stage 4: gene level ----------------
def gene_level(fam):
    from digas_raw import build
    D = build("exon_5k"); iids = D["iids"]
    grp = fam.set_index("iid").group.reindex(iids).values
    V = region_values(D, np.ones(len(iids), bool))                       # genes x subjects (cohort minor allele)
    nsnp = np.asarray((D["Wc"] > 0).sum(0)).ravel()
    rng = np.random.default_rng(SEED)
    H, B = {}, {}
    for g in GROUPS:
        i = np.where(grp == g)[0]; Vi = V[:, i].T                         # subjects x genes
        H[g] = Vi.mean(0)
        w = rng.multinomial(len(i), np.full(len(i), 1 / len(i)), size=NBOOT).astype(np.float32) / len(i)
        B[g] = w @ Vi                                                    # NBOOT x genes
    G = pd.DataFrame(dict(gene=D["genes"], n_snps=nsnp, **{f"hill_{g}": H[g] for g in GROUPS}))
    G["label"] = G.gene.str.extract(r"(ENSG\d+)")[0].map(KNOWN).fillna("")
    for a, b in (("CN", "AD"), ("CN", "MCI")):
        d = np.log1p(H[b]) - np.log1p(H[a]); db = np.log1p(B[b]) - np.log1p(B[a])
        se = db.std(0, ddof=1); k = f"{a}_vs_{b}"
        G[f"FC_{k}"] = np.abs(d); G[f"se_{k}"] = se
        with np.errstate(divide="ignore", invalid="ignore"): G[f"z_{k}"] = np.where(se > 0, d / se, np.nan)
        G[f"rank_FC_{k}"] = G[f"FC_{k}"].rank(ascending=False, method="min")
    rows = []
    bins = [0, 1, 2, 5, 10, 20, 50, 10**6]
    for k in ("CN_vs_AD", "CN_vs_MCI"):
        G["nbin"] = pd.cut(G.n_snps, bins, labels=["1", "2", "3-5", "6-10", "11-20", "21-50", ">50"])
        top = G.nsmallest(1270, f"rank_FC_{k}")
        for nb_, g in G.groupby("nbin", observed=True):
            z = g[f"z_{k}"].dropna()
            rows.append(dict(comparison=k.replace("_vs_", " vs "), snps_per_gene=nb_, genes=len(g),
                             share_of_genes=len(g) / len(G), share_of_top1270_FC=(top.nbin == nb_).mean(),
                             FC_median=g[f"FC_{k}"].median(), SE_median=g[f"se_{k}"].median(),
                             z_sd=z.std(), share_abs_z_gt_2=(z.abs() > 2).mean(), expected_abs_z_gt_2=0.0455))
    S = pd.DataFrame(rows)
    G.drop(columns="nbin").to_csv(OUT / "gene_variability_exon_5k.csv.gz", index=False)
    S.to_csv(OUT / "gene_variability_by_snp_count.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    k = "CN_vs_AD"; x = G.n_snps.clip(upper=200)
    ax[0].hexbin(x, G[f"se_{k}"], xscale="log", yscale="log", gridsize=45, bins="log", cmap="Blues", mincnt=1)
    ax[0].set_xlabel("SNPs in gene (exon_5k, log)"); ax[0].set_ylabel("bootstrap SE of CN-AD log Hill difference")
    ax[0].set_title("Sampling noise of the DiGAS FC shrinks with SNP count", loc="left", fontsize=9, color=INK)
    z = G[f"z_{k}"].dropna().clip(-8, 8)
    ax[1].hist(z, bins=120, density=True, color=COL["CN"], alpha=.75, label="observed genes")
    xx = np.linspace(-6, 6, 300); ax[1].plot(xx, stats.norm.pdf(xx), color=INK, lw=1.5, label="N(0,1): sampling only")
    for _, r in G[G.label != ""].iterrows():
        if np.isfinite(r[f"z_{k}"]):
            ax[1].axvline(r[f"z_{k}"], color=COL["AD"], lw=1); ax[1].text(r[f"z_{k}"], .38, r.label.split(" ")[0],
                                                                       rotation=90, fontsize=7, color=MUTED, va="top")
    ax[1].set_xlabel("z = (AD - CN) / bootstrap SE, per gene"); ax[1].set_ylabel("density")
    ax[1].set_title("CN vs AD gene differences vs sampling noise", loc="left", fontsize=9, color=INK)
    ax[1].legend(frameon=False, fontsize=8, loc="upper left")
    savefig(fig, "gene_noise_exon_5k")
    return S, G


def plots_snp(D, pruned, P):
    Dp = D[D.snp.isin(pruned)]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    for g in GROUPS:                                             # SE vs MAF
        x = np.linspace(.01, .5, 100); n = (D[f"n_{g}"]).median()
        ax[0].plot(x, np.sqrt(x * (1 - x) / (2 * n)), color=COL[g], lw=2, label=f"{g} (n={int(n)})")
    ax[0].set_xlabel("minor allele frequency"); ax[0].set_ylabel("standard error of the group frequency")
    ax[0].set_title("Sampling error of each group's allele frequency", loc="left", fontsize=9, color=INK)
    ax[0].legend(frameon=False, fontsize=8)
    for k, (a, b) in enumerate([("CN", "AD"), ("CN", "AD+MCI")], start=1):   # QQ, pruned set
        p = np.sort(Dp[f"p_{a}_vs_{b}"].dropna().values); e = (np.arange(1, len(p) + 1) - .5) / len(p)
        ax[k].plot(-np.log10(e), -np.log10(p), ".", ms=2.5, color=COL["CN"], label="LD-pruned SNPs")
        p2 = np.sort(D[f"p_{a}_vs_{b}"].dropna().values); e2 = (np.arange(1, len(p2) + 1) - .5) / len(p2)
        ax[k].plot(-np.log10(e2), -np.log10(p2), ".", ms=2, color=COL["MCI"], alpha=.5, label="all SNPs")
        m = -np.log10(e2[0]); ax[k].plot([0, m], [0, m], color=INK, lw=1, ls="--", label="sampling only")
        ax[k].set_xlabel("expected -log10 p (no true difference)"); ax[k].set_ylabel("observed -log10 p")
        ax[k].set_title(f"{a} vs {b}: allele frequency differences", loc="left", fontsize=9, color=INK)
        ax[k].legend(frameon=False, fontsize=8, loc="upper left")
    savefig(fig, "snp_noise_qq")

    fig, ax = plt.subplots(figsize=(6.5, 3.9))
    for (a, b), c in zip(PAIRS, ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]):
        q = P[(P.comparison == f"{a} vs {b}") & (P.alpha == 0.05)]
        ax.plot(range(len(q)), q.min_detectable_OR_median, "-o", color=c, lw=2, ms=5, label=f"{a} vs {b}, alpha 0.05")
        q = P[(P.comparison == f"{a} vs {b}") & (P.alpha == 5e-8)]
        ax.plot(range(len(q)), q.min_detectable_OR_median, "--o", color=c, lw=1.5, ms=4, alpha=.8)
    ax.axhspan(1.1, 1.25, color="#ecebe7", zorder=0); ax.text(0, 1.17, "typical AD risk loci (OR 1.1-1.25)",
                                                             fontsize=7, color=MUTED, va="center")
    ax.set_xticks(range(len(q)), q.maf_bin); ax.set_xlabel("minor allele frequency"); ax.set_yscale("log")
    ax.set_yticks([1.1, 1.25, 1.5, 2, 3, 5]); ax.set_yticklabels(["1.1", "1.25", "1.5", "2", "3", "5"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylabel("smallest detectable OR (80% power)")
    ax.set_title("Effect sizes this dataset can see (solid: alpha 0.05, dashed: 5e-8)", loc="left", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=7)
    savefig(fig, "power_min_detectable_or")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fam, bim = load_meta()
    print("group sizes:", fam.group.value_counts().to_dict(), flush=True)
    print("stage 0: structure", flush=True); S0, kin = structure(fam)
    pruned = set(pd.read_csv(TMP / "prune.prune.in", header=None)[0])
    print("stage 1-2: per SNP", flush=True); D = per_snp(fam, bim)
    S2 = cohort_summary(D, pruned); S2.to_csv(OUT / "cohort_summary.csv", index=False)
    print("stage 3: power", flush=True); P, need = power_tables(D, fam)
    P.to_csv(OUT / "power_min_detectable_or.csv", index=False); need.to_csv(OUT / "power_sample_size_needed.csv", index=False)
    top = D.nsmallest(30, "p_CN_vs_AD+MCI")[["snp", "chr", "bp"] + [f"maf_{g}" for g in GROUPS] +
                                           ["p_CN_vs_AD", "p_CN_vs_AD+MCI", "p_3groups", "minOR_CN_vs_AD+MCI_5e-08"]].copy()
    top["genes_within_20kb"] = gene_map(top); top.to_csv(OUT / "top_snps_CN_vs_ADMCI.csv", index=False)
    D.to_csv(OUT / "snp_by_group.csv.gz", index=False)
    plots_snp(D, pruned, P)
    print("stage 4: gene level", flush=True); S4, G = gene_level(fam)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    for name, T in (("structure", S0), ("kinship", kin), ("cohort", S2), ("power", need), ("genes", S4)):
        print(f"\n== {name}\n{T.round(4).to_string(index=False)}")
    print("\n== known genes\n", G[G.label != ""][["label", "n_snps", "FC_CN_vs_AD", "z_CN_vs_AD", "rank_FC_CN_vs_AD",
                                                  "z_CN_vs_MCI"]].round(3).to_string(index=False))
    print("\n== top SNPs\n", top.head(15).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
