#!/usr/bin/env python3
"""Which regions do DiGAS and SKAT select in common? Qualitative part (needs selection_overlap.py outputs).

Annotation (downloaded to annotation/, public files):
  gencode.v44lift37.basic.annotation.gtf.gz   gene symbol, biotype (GRCh37 coordinates; matched by stable ENSG ID)
  hgnc_complete_set.txt                       gene description
  goa_human.gaf.gz + go-basic.obo             GO annotations, propagated up is_a / part_of
  ReactomePathways.gmt.zip                    Reactome pathways (symbols)
  gwas_ad/page_*.json                         GWAS Catalog v2 associations for Alzheimer disease (MONDO_0004975)

Per case (exon_5k ADNI1/2/3, ADNI1 exon / exon_protein_coding) and set (Common, DiGAS-only, SKAT-only):
  1. biotype composition vs all tested regions
  2. known AD genes (GWAS Catalog mapped genes, p < 5e-8): which are tested / selected, Fisher enrichment
  3. GO (BP/MF/CC) and Reactome over-representation, universe = tested regions with a symbol; one-sided Fisher,
     BH FDR within each library; terms with 10-500 tested genes
  4. a per-region table of the common set with symbol, biotype, description, p-values
Caveat: exon_5k regions overlap and neighbouring genes share SNPs, so clustered gene families can inflate enrichment;
we report the number of distinct loci (genes > 1 Mb apart) behind each enriched term.

  .venv/bin/python scripts/analysis/selection_overlap_annot.py
Outputs: results/selection_overlap/annot/
"""
import glob, gzip, io, json, re, zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
import importlib.util
_s = importlib.util.spec_from_file_location("O", ROOT / "scripts/analysis/selection_overlap.py")
O = importlib.util.module_from_spec(_s); _s.loader.exec_module(O)
A = ROOT / "annotation"; IN = ROOT / "results/selection_overlap"; OUT = IN / "annot"; OUT.mkdir(parents=True, exist_ok=True)
CASES = [("ADNI1", "exon_5k"), ("ADNI2", "exon_5k"), ("ADNI3", "exon_5k"), ("ADNI1", "exon"), ("ADNI1", "exon_protein_coding")]
stable = lambda g: g.split(".")[0]


def gene_table():
    rows = []
    with gzip.open(A / "gencode.v44lift37.basic.annotation.gtf.gz", "rt") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.split("\t", 8)
            if p[2] != "gene":
                continue
            at = dict(re.findall(r'(\S+) "([^"]*)"', p[8]))
            rows.append((stable(at["gene_id"]), at.get("gene_name"), at.get("gene_type"), p[0].replace("chr", ""),
                         int(p[3]), int(p[4])))
    g = pd.DataFrame(rows, columns=["ENSG", "Symbol", "Biotype", "chr", "start", "end"]).drop_duplicates("ENSG").set_index("ENSG")
    h = pd.read_csv(A / "hgnc_complete_set.txt", sep="\t", usecols=["ensembl_gene_id", "name", "symbol"], dtype=str)
    h = h.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
    g["Description"] = h.name.reindex(g.index)
    miss = h.index.difference(g.index)          # HGNC genes absent from GENCODE v44 (retired IDs): keep symbol
    g = pd.concat([g, pd.DataFrame({"Symbol": h.loc[miss, "symbol"], "Description": h.loc[miss, "name"]})])
    return g


def go_sets(symbols):
    parents = defaultdict(set); names, ns = {}, {}
    cur = None
    for ln in open(A / "go-basic.obo"):
        ln = ln.strip()
        if ln == "[Term]":
            cur = {}
        elif ln.startswith("[") and ln != "[Term]":
            cur = None
        elif cur is not None and ":" in ln:
            k, v = ln.split(": ", 1)
            if k == "id":
                cur["id"] = v
            elif k == "name":
                names[cur["id"]] = v
            elif k == "namespace":
                ns[cur["id"]] = {"biological_process": "GO_BP", "molecular_function": "GO_MF", "cellular_component": "GO_CC"}[v]
            elif k == "is_a":
                parents[cur["id"]].add(v.split()[0])
            elif k == "relationship" and v.startswith("part_of"):
                parents[cur["id"]].add(v.split()[1])
    anc = {}

    def ancestors(t):
        if t not in anc:
            s = {t}
            for p in parents.get(t, ()):
                s |= ancestors(p)
            anc[t] = s
        return anc[t]

    import sys; sys.setrecursionlimit(10000)
    sets = defaultdict(set)
    with gzip.open(A / "goa_human.gaf.gz", "rt") as f:
        for ln in f:
            if ln.startswith("!"):
                continue
            p = ln.split("\t")
            if "NOT" in p[3] or p[2] not in symbols:
                continue
            for t in ancestors(p[4]):
                sets[t].add(p[2])
    out = defaultdict(dict)
    for t, s in sets.items():
        if t in ns:
            out[ns[t]][f"{t} {names.get(t, '')}"] = s
    return out


def reactome_sets(symbols):
    z = zipfile.ZipFile(A / "ReactomePathways.gmt.zip")
    sets = {}
    for ln in io.TextIOWrapper(z.open(z.namelist()[0]), encoding="utf-8"):
        p = ln.rstrip("\n").split("\t")
        sets[f"{p[1]} {p[0]}"] = set(p[2:]) & symbols
    return sets


def gwas_ad_genes(min_papers=3):
    """GWAS Catalog mapped genes of AD associations with p < 5e-8. The AD trait also pulls in proxy/family-history,
    pleiotropy and PRS studies (1,969 genes, 1,700 from a single paper), so 'known AD genes' = genes reported in
    >= min_papers independent publications (135 genes: APOE, BIN1, CLU, PICALM, CR1, TREM2, SORL1, ...)."""
    genes = defaultdict(lambda: {"min_p": 1.0, "n_assoc": 0, "papers": set()})
    for f in glob.glob(str(A / "gwas_ad/page_*.json")):
        for a in json.load(open(f)).get("_embedded", {}).get("associations", []):
            p = a.get("p_value") or 1.0
            if p >= 5e-8:
                continue
            for g in a.get("mapped_genes") or []:
                for s in re.split(r"[ ,;]+", g):
                    if s and s not in ("-", "Y_RNA"):
                        d = genes[s]; d["min_p"] = min(d["min_p"], p); d["n_assoc"] += 1; d["papers"].add(a.get("pubmed_id"))
    t = pd.DataFrame({s: dict(min_p=d["min_p"], n_assoc=d["n_assoc"], n_papers=len(d["papers"])) for s, d in genes.items()}).T
    t["Replicated"] = t.n_papers >= min_papers
    return t


def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p)
    q = p[o] * m / np.arange(1, m + 1); q = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    out = np.empty(m); out[o] = q
    return out


def mh_greater(hits, s, strata):
    """size-adjusted one-sided test: under the null each term gene is selected with the selection rate of its size
    class (SNP-count bin), independently; P(X >= observed hits) from the exact Poisson-binomial distribution."""
    probs = []
    for members in strata.values():
        n = len(members)
        if n == 0:
            continue
        r = len(members & hits) / n
        probs += [r] * len(members & s)
    k = len(s & hits & set().union(*strata.values()))
    dist = np.zeros(len(probs) + 1); dist[0] = 1.0
    for p in probs:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p; dist[0] *= (1 - p)
    return float(dist[k:].sum()), float(np.sum(probs))


def enrich(hits, universe, libs, strata, min_size=10, max_size=500, gpos=None):
    res = []
    for lib, sets in libs.items():
        rows = []
        for term, s in sets.items():
            s = s & universe
            if not (min_size <= len(s) <= max_size):
                continue
            k = len(s & hits)
            if k < 2:
                continue
            p = fisher_exact([[k, len(hits) - k], [len(s) - k, len(universe) - len(hits) - len(s) + k]], alternative="greater")[1]
            pa, ea = mh_greater(hits, s, strata)
            rows.append(dict(Library=lib, Term=term, Size=len(s), Hits=k, Expected=len(s) * len(hits) / len(universe),
                             p=p, Expected_size_adj=ea, p_size_adj=pa, Genes=",".join(sorted(s & hits)),
                             Loci=n_loci(s & hits, gpos)))
        d = pd.DataFrame(rows)
        if len(d):
            d["FDR"] = bh(d.p); d["FDR_size_adj"] = bh(d.p_size_adj)
            res.append(d.sort_values("p_size_adj"))
    return pd.concat(res) if res else pd.DataFrame()


def n_loci(genes, gpos):
    if gpos is None or not genes:
        return np.nan
    p = gpos.loc[gpos.Symbol.isin(genes)].dropna(subset=["start"]).sort_values(["chr", "start"])
    n, last = 0, (None, -1e12)
    for c, s in zip(p.chr, p.start):
        if c != last[0] or s - last[1] > 1e6:
            n += 1
        last = (c, s)
    return n


def main():
    G = gene_table()
    ad = gwas_ad_genes(); ad.to_csv(OUT / "gwas_ad_genes.csv")
    syms_all = set(G.Symbol.dropna())
    libs = {**go_sets(syms_all), "Reactome": reactome_sets(syms_all)}
    gpos = G.dropna(subset=["Symbol"]).drop_duplicates("Symbol")[["Symbol", "chr", "start"]]
    bt_rows, ad_rows, en_all, summ, agree_rows, sexchr_rows = [], [], [], [], [], []
    # sex x diagnosis per cohort
    sexdx = []
    for ds in ("ADNI1", "ADNI2", "ADNI3"):
        fam = pd.read_csv(ROOT / f"{ds}_QC.fam", sep=r"\s+", header=None, usecols=[1, 4], names=["IID", "Sex"], dtype=str)
        ph = pd.read_csv(ROOT / f"pheno/{ds.lower()}pheno.txt", sep=r"\s+", header=None, usecols=[0, 1], names=["IID", "Dx"], dtype=str)
        m = fam.merge(ph, on="IID")
        for dx in ("CN", "MCI", "AD"):
            g = m[m.Dx == dx]; k = g[g.Sex.isin(["1", "2"])]
            sexdx.append(dict(Dataset=ds, Dx=dx, N=len(g), Sex_unknown=int((~g.Sex.isin(["1", "2"])).sum()),
                              Male_frac=(k.Sex == "1").mean()))
        mk = m[m.Sex.isin(["1", "2"])]
        for other in ("AD", "MCI"):
            t = pd.crosstab(mk[mk.Dx.isin(["CN", other])].Dx, mk[mk.Dx.isin(["CN", other])].Sex)
            sexdx.append(dict(Dataset=ds, Dx=f"CN vs {other} chi2 p", N=np.nan, Male_frac=chi2_contingency(t)[1]))
    sexdx = pd.DataFrame(sexdx); sexdx.to_csv(OUT / "sex_by_diagnosis.csv", index=False)
    for ds, reg in CASES:
        sel = pd.read_csv(IN / f"selected_{ds}_{reg}.csv")
        tested = pd.read_csv(IN / f"common_{ds}_{reg}.csv")   # common set table
        # universe = regions tested by both tools in this cohort (same as selection_overlap.py)
        pdg, psk = O.CASES[(ds, reg)]
        Sall = O.load_all(psk)
        uni_ids = np.array(sorted(O.load_all(pdg).index.intersection(Sall.index)))
        nsnp = Sall.loc[uni_ids, "N_snps"]
        ann = lambda ids: G.reindex([stable(i) for i in ids])
        D = set(sel.Gene[sel.DiGAS]); S = set(sel.Gene[sel.SKAT])
        groups = {"Common": D & S, "DiGAS-only": D - S, "SKAT-only": S - D}
        # agreement by region size
        for lab, g in nsnp.groupby(pd.cut(nsnp, [0, 1, 2, 5, 10, 1e9], labels=["1", "2", "3-5", "6-10", ">10"])):
            ids = set(g.index); a_, b_ = D & ids, S & ids
            agree_rows.append(dict(Dataset=ds, Region=reg, SNPs=lab, Regions=len(ids), DiGAS=len(a_), SKAT=len(b_),
                                   Common=len(a_ & b_), Common_of_DiGAS=len(a_ & b_) / max(len(a_), 1),
                                   Common_of_SKAT=len(a_ & b_) / max(len(b_), 1)))
        for gname, ids in [("All tested", uni_ids)] + list(groups.items()):
            v = nsnp.reindex(list(ids))
            agree_rows.append(dict(Dataset=ds, Region=reg, SNPs=f"[{gname}] median", Regions=len(ids),
                                   Common_of_DiGAS=np.nan, DiGAS=v.median(), SKAT=(v == 1).mean()))
        U_ = ann(uni_ids)
        summ.append(dict(Dataset=ds, Region=reg, Regions=len(uni_ids), Annotated=int(U_.Symbol.notna().sum())))
        # 1. biotypes
        top_bt = U_.Biotype.value_counts().index[:6]
        b = lambda x: x.where(x.isin(top_bt), "other").fillna("unannotated")
        for gname, ids in [("All tested", uni_ids)] + list(groups.items()):
            vc = b(ann(ids).Biotype).value_counts(normalize=True)
            bt_rows.append(dict(Dataset=ds, Region=reg, Set=gname, N=len(ids), **vc.round(4).to_dict()))
        # 1b. sex chromosomes (neither tool adjusts for sex)
        chrom = pd.Series(U_.chr.astype(str).values, index=uni_ids)     # GENCODE chromosome (region files code X and Y as 23)
        for gname, ids in groups.items():
            for cx in ("X", "Y"):
                inset = chrom.loc[sorted(ids)] == cx; allc = chrom == cx
                k, n = int(inset.sum()), len(ids); K, N = int(allc.sum()), len(chrom)
                p = fisher_exact([[k, n - k], [K - k, N - n - K + k]], alternative="greater")[1]
                sexchr_rows.append(dict(Dataset=ds, Region=reg, Set=gname, Chr=cx, In_set=k,
                                        Set_size=n, Frac_set=k / n, Frac_all=K / N, Fold=(k / n) / (K / N) if K else np.nan,
                                        p=p))
        # 2. known AD genes
        sym_of = dict(zip(uni_ids, U_.Symbol.values))
        uni_sym = {s for s in sym_of.values() if isinstance(s, str)}
        ad_t = set(ad.index[ad.Replicated]) & uni_sym
        # gene size class = SNP-count bin of its largest region (used by the size-adjusted tests)
        sz = pd.Series(nsnp.values, index=[sym_of[i] for i in uni_ids]).dropna()
        sz = sz[[isinstance(i, str) for i in sz.index]].groupby(level=0).max()
        sbin = pd.cut(sz, [0, 1, 2, 3, 5, 10, 20, 1e9], labels=False)
        strata = {b_: set(g.index) for b_, g in sz.groupby(sbin)}
        for gname, ids in groups.items():
            hs = {sym_of[i] for i in ids if isinstance(sym_of.get(i), str)}
            k = len(hs & ad_t)
            p = fisher_exact([[k, len(hs) - k], [len(ad_t) - k, len(uni_sym) - len(hs) - len(ad_t) + k]], alternative="greater")[1]
            pa, ea = mh_greater(hs, ad_t, strata)
            ad_rows.append(dict(Dataset=ds, Region=reg, Set=gname, Genes=len(hs), AD_genes_tested=len(ad_t), AD_hits=k,
                                Loci=n_loci(hs & ad_t, gpos), AD_loci_tested=n_loci(ad_t, gpos),
                                Expected=len(hs) * len(ad_t) / len(uni_sym), p=p, Expected_size_adj=ea, p_size_adj=pa,
                                Hits=",".join(sorted(hs & ad_t))))
        # 3. enrichment
        for gname, ids in groups.items():
            hs = {sym_of[i] for i in ids if isinstance(sym_of.get(i), str)}
            e = enrich(hs, uni_sym, libs, strata, gpos=gpos)
            if len(e):
                e.insert(0, "Set", gname); e.insert(0, "Region", reg); e.insert(0, "Dataset", ds); en_all.append(e)
        # 4. common-set table
        t = tested.copy(); a = ann(t.Gene)
        for c in ("Symbol", "Biotype", "chr", "start", "end", "Description"):
            t[c] = a[c].values
        t["Known_AD_GWAS"] = t.Symbol.isin(ad.index[ad.Replicated])
        t["AD_GWAS_any"] = t.Symbol.isin(ad.index)
        t["best_p"] = t.filter(like="_p_").min(axis=1)
        t.sort_values("best_p").to_csv(OUT / f"common_annotated_{ds}_{reg}.csv", index=False)

    ag = pd.DataFrame(agree_rows)
    ag[~ag.SNPs.str.startswith("[")].to_csv(IN / "agreement_by_size.csv", index=False)
    ag[ag.SNPs.str.startswith("[")].rename(columns={"DiGAS": "Median_SNPs", "SKAT": "Frac_1SNP"})[
        ["Dataset", "Region", "SNPs", "Regions", "Median_SNPs", "Frac_1SNP"]].to_csv(IN / "set_sizes.csv", index=False)
    sx = pd.DataFrame(sexchr_rows); sx.to_csv(OUT / "sex_chromosomes.csv", index=False)
    bt = pd.DataFrame(bt_rows).fillna(0); bt.to_csv(OUT / "biotypes.csv", index=False)
    adr = pd.DataFrame(ad_rows); adr.to_csv(OUT / "known_ad_genes.csv", index=False)
    en = pd.concat(en_all) if en_all else pd.DataFrame(); en.to_csv(OUT / "enrichment.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    rep = ["== annotation coverage", pd.DataFrame(summ).to_string(index=False), "",
           f"== GWAS Catalog AD genes (p < 5e-8): {len(ad)} mapped; {int(ad.Replicated.sum())} in >= 3 publications "
           "(used as known AD genes)", "",
           "== agreement by region size", ag.round(3).to_string(index=False), "",
           "== sex by diagnosis", sexdx.round(4).to_string(index=False), "",
           "== sex-chromosome regions in each set", sx.round(4).to_string(index=False), "",
           "== biotypes", bt.to_string(index=False), "", "== known AD genes", adr.round(4).to_string(index=False), ""]
    if len(en):
        cols = ["Dataset", "Region", "Set", "Library", "Term", "Size", "Hits", "Expected", "Expected_size_adj", "Loci",
                "FDR", "FDR_size_adj"]
        sig = en[(en.FDR < 0.05) | (en.FDR_size_adj < 0.05)]
        rep += [f"== enrichment: terms with FDR < 0.05 (unadjusted: {(en.FDR < .05).sum()}, size-adjusted: "
                f"{(en.FDR_size_adj < .05).sum()})", sig[cols].round(4).to_string(index=False) if len(sig) else "(none)", ""]
        # replication of terms across cohorts (exon_5k), nominal p_size_adj < 0.01
        e5 = en[(en.Region == "exon_5k") & (en.p_size_adj < 0.01)]
        rep_t = e5.groupby(["Set", "Library", "Term"]).Dataset.nunique()
        rep += ["== exon_5k terms with size-adjusted p < 0.01 in 2+ cohorts (same set)",
                rep_t[rep_t >= 2].to_string() if (rep_t >= 2).any() else "(none)", "",
                "== best 3 terms per case/set (size-adjusted)",
                en.groupby(["Dataset", "Region", "Set"]).head(3)[["Dataset", "Region", "Set", "Library", "Term", "Hits",
                                                                   "Expected_size_adj", "Loci", "p_size_adj",
                                                                   "FDR_size_adj"]].round(4).to_string(index=False)]
    open(OUT / "report.txt", "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
