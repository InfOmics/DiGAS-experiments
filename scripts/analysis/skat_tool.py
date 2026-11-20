#!/usr/bin/env python3
"""SKAT drop-in replacement for scripts/nested/digas3-linux_validation.py (same command line, same output files).

Lets the existing DiGAS pipelines run unchanged with SKAT as the region-selection step:
  nested:  scripts/nested/digas_classification_nested_10fold.py --digas-script scripts/analysis/skat_tool.py ...
  leaky:   run this tool on all subjects, then the nested pipeline with scripts/analysis/precomputed_selection.py
           (see scripts/analysis/tool_comparison_queue.py)

SKAT settings, as used in the paper:
  SKAT (not SKAT-O), default weights Beta(1,25), linear-weighted kernel, no covariates, all regions of the region file,
  binary phenotype per contrast (CN = control). Contrasts CN vs AD and CN vs MCI, as for DiGAS (only those whose groups
  are present in --pheno). R package SKAT defaults otherwise: SKAT_Null_Model(y ~ 1, out_type = "D"), which applies the
  small-sample adjustment for n < 2000; missing genotypes imputed by the mean ("fixed"); alleles flipped to the minor one.
Region -> SNP assignment: same strict interval rule as DiGAS and the nested script (start < bp < end, same chromosome);
a gene spans several rows (one per exon); each SNP enters its gene's set once.

Outputs (same names/columns as DiGAS): <out>_all.txt and <out>_significant.txt (p <= 0.05), space separated,
columns Cat1 Cat2 Gene FC Pvalue HillCat1 HillCat2 (FC/Hill = NA; extra columns N_markers N_markers_tested).
Usage: python skat_tool.py --bfile ADNI1_QC --region region/exon_5k.txt --pheno pheno.txt --out prefix [--iteration N]
       (--iteration is accepted and ignored)
"""
import argparse, os, subprocess, tempfile, shutil
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RSCRIPT = ROOT / ".r-env" / "R" / "bin" / "Rscript"
PLINK = ROOT / "plink"
CACHE = ROOT / "results" / "skat_setid"
CONTRASTS = [("CN", "AD"), ("CN", "MCI")]

R_CODE = r"""
suppressMessages(library(SKAT))
a <- commandArgs(trailingOnly = TRUE); pre <- a[1]; setid <- a[2]; out <- a[3]
ssd <- paste0(pre, ".SSD"); info <- paste0(pre, ".SSD.info")
invisible(capture.output(Generate_SSD_SetID(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"),
                                            setid, ssd, info)))
FAM <- Read_Plink_FAM(paste0(pre, ".fam"), Is.binary = TRUE)
y <- FAM$Phenotype
set.seed(1)  # small-sample adjustment (n < 2000) resamples; fixed seed -> reproducible p-values
obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
SSD.INFO <- Open_SSD(ssd, info)
res <- SKAT.SSD.All(SSD.INFO, obj)$results
Close_SSD()
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
"""


def setid_file(region_file, bim_file):
    """SetID file (region<TAB>SNP) for the region file, cached per (dataset, region file).
    The SNP content depends on the .bim, so each dataset has its own file (ADNI1 keeps the old unprefixed name)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    ds = Path(bim_file).stem
    f = CACHE / (f"{Path(region_file).stem}.SetID" if ds == "ADNI1_QC" else f"{ds}_{Path(region_file).stem}.SetID")
    if f.exists():
        return f
    bim = pd.read_csv(bim_file, sep=r"\s+", header=None, usecols=[0, 1, 3], names=["chr", "snp", "bp"],
                      dtype={"chr": str, "snp": str})
    rows = []
    for line in open(region_file):
        if not line.strip() or line.startswith("#"):
            continue
        x = line.rstrip().split("\t")
        if len(x) == 2:
            x = [x[0] + "-" + x[1], x[0], x[1]]
        c1, s = x[1].split(":", 1); c2, e = x[2].split(":", 1); assert c1 == c2
        rows.append((x[0], str(int(float(c1))) if c1.isdigit() else c1, int(s), int(e)))
    out = []
    for c, g in pd.DataFrame(rows, columns=["gene", "chr", "s", "e"]).groupby("chr", sort=False):
        b = bim[bim.chr == c].sort_values("bp", kind="stable"); bp = b.bp.values; ids = b.snp.values
        lo = np.searchsorted(bp, g.s.values, side="right"); hi = np.searchsorted(bp, g.e.values, side="left")
        for gene, l, h in zip(g.gene.values, lo, hi):
            for snp in dict.fromkeys(ids[l:h]):                     # dedupe within region, keep order
                out.append((gene, snp))
    S = pd.DataFrame(out, columns=["gene", "snp"]).drop_duplicates().sort_values("gene", kind="stable")
    tmp = f.with_suffix(".tmp"); S.to_csv(tmp, sep="\t", header=False, index=False); tmp.rename(f)
    return f


def run_contrast(bfile, pheno, ctrl, case, setid, work):
    sub = pheno[pheno.group.isin([ctrl, case])]
    fam = pd.read_csv(f"{bfile}.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["fid", "iid"], dtype=str)
    sub = fam.merge(sub, on="iid")
    sub["y"] = np.where(sub.group == case, 2, 1)
    sub[["fid", "iid"]].to_csv(work / "keep.txt", sep=" ", header=False, index=False)
    sub[["fid", "iid", "y"]].to_csv(work / "pheno.txt", sep=" ", header=False, index=False)
    pre = work / f"{ctrl}_{case}"
    subprocess.run([str(PLINK), "--bfile", bfile, "--keep", str(work / "keep.txt"), "--pheno", str(work / "pheno.txt"),
                    "--allow-no-sex", "--make-bed", "--out", str(pre)], check=True, stdout=subprocess.DEVNULL)
    (work / "skat.R").write_text(R_CODE)
    res = work / f"{ctrl}_{case}.tsv"
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")    # 1 thread is ~2x faster than 8 here
    subprocess.run([str(RSCRIPT), str(work / "skat.R"), str(pre), str(setid), str(res)], check=True,
                   stdout=subprocess.DEVNULL, env=env)
    r = pd.read_csv(res, sep="\t", dtype={"SetID": str})
    return pd.DataFrame({"Cat1": ctrl, "Cat2": case, "Gene": r.SetID, "FC": "NA", "Pvalue": r["P.value"],
                         "HillCat1": "NA", "HillCat2": "NA", "N_markers": r["N.Marker.All"],
                         "N_markers_tested": r["N.Marker.Test"]})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bfile", required=True); p.add_argument("--region", required=True)
    p.add_argument("--pheno", required=True); p.add_argument("--out", required=True)
    p.add_argument("--iteration", type=int, default=0)
    p.add_argument("--saveTemporary", action="store_false")
    a = p.parse_args()
    pheno = pd.read_csv(a.pheno, sep=r"\s+", header=None, usecols=[0, 1], names=["iid", "group"], dtype=str)
    setid = setid_file(a.region, f"{a.bfile}.bim")
    work = Path(tempfile.mkdtemp(prefix="skat_", dir=str(Path(a.out).parent)))
    try:
        parts = [run_contrast(a.bfile, pheno, c, k, setid, work) for c, k in CONTRASTS
                 if (pheno.group == c).any() and (pheno.group == k).any()]
    finally:
        shutil.rmtree(work, ignore_errors=True)
    R = pd.concat(parts, ignore_index=True)
    R.to_csv(f"{a.out}_all.txt", sep=" ", index=False, na_rep="NA")
    R[pd.to_numeric(R.Pvalue, errors="coerce") <= 0.05].to_csv(f"{a.out}_significant.txt", sep=" ", index=False,
                                                                na_rep="NA")
    print({f"{c}_vs_{k}": int((pd.to_numeric(g.Pvalue) <= 0.05).sum()) for (c, k), g in R.groupby(["Cat1", "Cat2"])})


if __name__ == "__main__":
    main()
