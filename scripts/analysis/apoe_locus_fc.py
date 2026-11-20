"""DiGAS raw FC (no p-values) for the APOE locus genes: APOE, TOMM40, APOC1 (+ ENSG00000267282, which holds rs2075650).

A. Each nested-CV training fold (saved DiGAS output): FC and its rank among all 44,170 genes, CN vs AD and CN vs MCI.
B. Whole cohort: FC and its rank among all genes with the real labels (results/apoe_locus_fc/cohort_fc.csv).

The training folds are those of the nested DiGAS exon_5k run of scripts/analysis/tool_comparison_queue.py
(results/tool_comparison/nested_DiGAS_exon_5k).
Outputs: results/apoe_locus_fc/{fold_fc.csv, cohort_fc.csv}
Usage: python scripts/analysis/apoe_locus_fc.py
"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from digas_raw import build, hills, fc, ROOT

OUT = os.path.join(ROOT, "results", "apoe_locus_fc")
FOLDS = os.path.join(ROOT, "results/tool_comparison/nested_DiGAS_exon_5k/nested_work/ADNI1_QC_exon_5k")
GENES = {"ENSG00000130203": "APOE", "ENSG00000130204": "TOMM40", "ENSG00000130208": "APOC1",
         "ENSG00000267282": "ENSG00000267282 (rs2075650)"}
CONTRASTS = [("CN", "AD"), ("CN", "MCI")]


def short(g): return g.split(".")[0]


def rank_of(v, idx):
    """1 = largest FC. Ties share the best rank."""
    return np.array([(v > v[i]).sum() + 1 for i in idx])


def part_a():
    rows = []
    for f in sorted(glob.glob(os.path.join(FOLDS, "fold_*/digas_train_all.txt"))):
        d = pd.read_csv(f, sep=" "); fold = os.path.basename(os.path.dirname(f))
        for a, b in CONTRASTS:
            x = d[((d.Cat1 == a) & (d.Cat2 == b)) | ((d.Cat1 == b) & (d.Cat2 == a))].reset_index(drop=True)
            idx = [i for i, g in enumerate(x.Gene) if short(g) in GENES]
            for i, r in zip(idx, rank_of(x.FC.values, idx)):
                rows.append(dict(fold=fold, contrast=f"{a} vs {b}", gene=GENES[short(x.Gene[i])], FC=x.FC[i],
                                 FC_rank=int(r), n_genes=len(x), FC_percentile=100 * (1 - (r - 1) / len(x)),
                                 Pvalue=x.Pvalue[i]))
    return pd.DataFrame(rows)


def group_fc(D, lab):
    H = {g: hills(D, lab == g) for g in ("CN", "MCI", "AD")}
    return {f"{a} vs {b}": fc(H[a], H[b]) for a, b in CONTRASTS}



def cohort_rank(D, lab):
    idx = [i for i, g in enumerate(D["genes"]) if short(g) in GENES]
    nsnp = np.asarray(D["Wc"].sum(0)).ravel(); rows = []
    for c, v in group_fc(D, lab).items():
        for i, r in zip(idx, rank_of(v, idx)):
            rows.append(dict(contrast=c, gene=GENES[short(D["genes"][i])], snps=int(nsnp[i]), FC=v[i], FC_rank=int(r)))
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True); pd.set_option("display.width", 220)
    A = part_a(); A.to_csv(os.path.join(OUT, "fold_fc.csv"), index=False)
    print("A) FC rank in each nested training fold (rank 1 = highest FC of 44,170 genes)")
    print(A.groupby(["contrast", "gene"]).agg(FC_median=("FC", "median"), rank_median=("FC_rank", "median"),
                                              rank_best=("FC_rank", "min"), rank_worst=("FC_rank", "max"),
                                              p_min=("Pvalue", "min"), p_max=("Pvalue", "max")).round(4).to_string())
    D = build(); lab = pd.Series(D["iids"]).map(pd.read_csv(f"{ROOT}/pheno/adni1pheno.txt", sep=r"\s+", header=None,
                                                             index_col=0)[1]).values.astype(object)
    B = cohort_rank(D, lab); B.to_csv(os.path.join(OUT, "cohort_fc.csv"), index=False)
    print("\nB) whole cohort, real labels"); print(B.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
