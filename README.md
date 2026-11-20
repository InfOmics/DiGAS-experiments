# DiGAS-experiments

Code and results of the additional analyses run in 2026 for the DiGAS paper, at the request of the editor and an additional reviewer.

The analyses cover:
- **Validation:** region selection inside cross-validation (nested) versus on the whole cohort (as in the paper).
- **Comparison with SKAT:** run through the same pipeline.
- **Unsupervised analysis:** PCA, clustering and embeddings in three ADNI cohorts.
- **Agreement:** the overlap between the regions selected by DiGAS and by SKAT.
- **Statistical power:** what the ADNI cohorts can detect.

Runs still in progress will be added when they finish: the other three region definitions of the paper, the supervised comparison on ADNI2 and ADNI3.

## Results and how to reproduce them

Run each script from the repository root. `bash replication/run_all.sh` runs them all in order.

| Result | Run script | Main code | Results in this repository |
|---|---|---|---|
| 0. Genotype QC (PLINK 1.9) | `scripts/qc/qc_plink19.sh <in> <out>` | — | — |
| 1. Classification with whole-cohort selection | `replication/1_supervised_comparison.sh` | `scripts/analysis/tool_comparison_queue.py`, `scripts/nested/`, `scripts/analysis/skat_tool.py`, `compare_tools.py` | `results/tool_comparison/summary_all.csv` (`arm = leaky`), `summary_folds.csv` |
| 2. Classification with nested selection | same run | same | `summary_all.csv` (`arm = nested`), `summary_paired_nested.csv` (fold-by-fold Wilcoxon test) |
| 3. Unsupervised analysis, ADNI1/2/3 | `replication/2_unsupervised_analysis.sh` | `scripts/analysis/unsupervised_tools.py`, `adni23_selections.py` | `results/unsupervised_tools/` (report: `report/unsupervised_DiGAS_vs_SKAT.pdf`) |
| 4. Agreement between DiGAS and SKAT | `replication/3_selection_overlap.sh` | `scripts/analysis/selection_overlap*.py` | `results/selection_overlap/` (report: `report/DiGAS_SKAT_common_selections.pdf`) |
| 4. APOE locus (APOC1) | `replication/4_apoc1_locus.sh` | `scripts/analysis/apoe_locus_fc.py` | `results/apoe_locus_fc/fold_fc.csv`, `cohort_fc.csv` |
| 5. Data quality and statistical power, ADNI1 | `replication/5_power_analysis.sh` | `scripts/analysis/snp_variability.py` | `results/snp_variability/` (report: `REPORT.md`) |

**Selected regions:** the regions each tool selected on the whole cohort (p ≤ 0.05), with their test statistic and p-value, are in:
- `results/tool_comparison/selection_all/`: ADNI1 exon and exon_protein_coding; ADNI2 and ADNI3 exon_5k.
- `results/fold_pair_scatter_leaky/digas_all_subjects/`: DiGAS, ADNI1 exon_5k.
- `results/skat_leaky/`: SKAT, ADNI1 exon_5k.

The file paths match those used in the document sent to the editor.

## Settings

- **Subjects:** ADNI1, 704 subjects of European ancestry after QC (CN 197, MCI 339, AD 168). The target is CN vs AD/MCI, as in the paper. ADNI2 (736 subjects) and ADNI3 (302 subjects with phenotype) are used in the unsupervised analysis.
- **Regions:** `region/exon.txt`, `region/exon_protein_coding.txt` and `region/exon_5k.txt` (each exon plus the 5 kb upstream of it, following the strand), from GENCODE v36lift37.
  - Each row is one interval; the first column is the region identifier (Ensembl gene ID).
  - Intervals with the same identifier form one SNP set.
- **Region selection:**
  - DiGAS: `scripts/nested/digas3-linux_validation.py`, 1,000 permutations.
  - SKAT: `scripts/analysis/skat_tool.py`, as in the paper: SKAT with default Beta(1,25) weights, no covariates.
  - A region is selected when its CN vs AD or CN vs MCI test has p ≤ 0.05.
- **Features:** every SNP in the selected regions, coded 1 if the subject carries the variant and 0 if not.
- **Classifiers** (scikit-learn): linear SVM, polynomial SVM, LDA and decision tree with default parameters, and kNN with k chosen by inner cross-validation on the training set.
- **Validation:** 10-fold stratified cross-validation with the same folds for both tools. Accuracy, macro-F1 and balanced accuracy are computed on the pooled out-of-fold predictions.
  - *Whole-cohort* setting: each tool selects regions once on all subjects, as in the paper.
  - *Nested* setting: selection is repeated inside each training fold.

## Data

The genotype and phenotype data come from the Alzheimer's Disease Neuroimaging Initiative (ADNI). They are not included and must be requested from ADNI (adni.loni.usc.edu). For the same reason, no file in `results/` contains subject identifiers or subject-level values; only summary results are included.

Expected layout, at the repository root:

```text
ADNI1_QC.bed / .bim / .fam                  PLINK files after QC (scripts/qc/qc_plink19.sh)
ADNI2/ADNI2_QC.*, ADNI3/ADNI3_QC.*         plus links in the root: ln -s ADNI2/ADNI2_QC.* ADNI3/ADNI3_QC.* .
pheno/adni1pheno.txt, adni2pheno.txt, adni3pheno.txt
                                            two columns, no header: subject ID (as in the .fam) and CN / MCI / AD
plink                                       PLINK 1.9 Linux binary (results produced with v1.90b6.25)
```

Result 4 also uses public annotation files: GENCODE, HGNC, Gene Ontology, Reactome and the GWAS Catalog. `scripts/analysis/download_annotation.sh` downloads them into `annotation/`; `3_selection_overlap.sh` calls it automatically. The published results used the versions of 18 November 2026. Later versions can change the enrichment results slightly. The list of known Alzheimer's genes that was used is in `results/selection_overlap/annot/gwas_ad_genes.csv`.

## Requirements

- **Linux:** the job queue reads `/proc/meminfo`.
- **Python 3.8 or later** in `.venv/` at the root: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`. The original runs used Python 3.8.10 with numpy 1.24.4, pandas 2.0.3, scikit-learn 1.3.2, scipy 1.10.1 and umap-learn 0.5.7.
- **R with the SKAT package:** expected at `.r-env/R/bin/Rscript`, the path is set in `scripts/analysis/skat_tool.py`. The original runs used R 4.3.3 and SKAT 2.2.5 from conda-forge (`r-base`, `r-skat`).
- **Memory and time:** 16 GB is enough with at most 2 jobs at once, which is the queue's default. The supervised comparison takes several days.
