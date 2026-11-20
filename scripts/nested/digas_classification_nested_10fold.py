#!/usr/bin/env python3
"""Validazione nested della pipeline DiGAS + classificazione su ADNI."""

import argparse
import math
import os
import time
import shutil
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, matthews_corrcoef,
    precision_recall_fscore_support, confusion_matrix,
)
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier


DATASETS = {
    "ADNI1_QC": "pheno/adni1pheno.txt",
    "ADNI2_QC": "pheno/adni2pheno.txt",
    "ADNI3_QC": "pheno/adni3pheno.txt",
}

METHODS = ["SVM_LINEAR", "SVM_POLY", "LDA", "DT", "KNN"]
SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description="Validazione nested DiGAS + classificazione")
    p.add_argument("--root", default=".")
    p.add_argument("--datasets", nargs="+", choices=DATASETS, default=["ADNI1_QC"])
    p.add_argument("--region-file", required=True)
    p.add_argument("--outer-folds", type=int, default=10)
    p.add_argument("--sampling-70-30", action="store_true")
    p.add_argument("--sampling-90-10", action="store_true")
    p.add_argument("--sampling-repeats", type=int, default=10)
    p.add_argument("--ad-vs-cn", action="store_true")
    p.add_argument("--digas-iterations", type=int, default=1000)
    p.add_argument("--outer-jobs", type=int, default=1)
    p.add_argument("--results-dir", default="classification_results")
    p.add_argument("--plink", default="./plink")
    p.add_argument("--digas-script", default="scripts/nested/digas3-linux_validation.py")
    p.add_argument("--keep-work", action="store_true")
    return p.parse_args()


def run_command(cmd, log_file, cwd=None):
    with open(log_file, "a") as log:
        log.write("\n$ " + " ".join(map(str, cmd)) + "\n")
        proc = subprocess.run(
            [str(x) for x in cmd],
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=cwd,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"Comando fallito. Vedi {log_file}")


def load_data(root, dataset, pheno_file, region_file, ad_vs_cn=False):
    # Fenotipo e soggetti disponibili nel FAM
    pheno = pd.read_csv(
        pheno_file, sep="\t", header=None, names=["IID", "OriginalGroup"], dtype=str
    )
    if ad_vs_cn:
        pheno = pheno[pheno["OriginalGroup"].isin(["CN", "AD"])].copy()
        pheno["Group"] = pheno["OriginalGroup"]
    else:
        pheno["Group"] = pheno["OriginalGroup"].map({"CN": "CN", "MCI": "AD", "AD": "AD"})

    fam = pd.read_csv(
        root / f"{dataset}.fam",
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
        dtype=str,
    )
    subjects = fam.merge(pheno, on="IID", how="inner").dropna().reset_index(drop=True)

    # BIM
    bim = pd.read_csv(
        root / f"{dataset}.bim",
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 3],
        names=["chr", "snp", "bp"],
        dtype={"chr": str, "snp": str},
    )

    def norm_chr(x):
        x = str(x).replace("chr", "").strip()
        aliases = {"X": "23", "Y": "24", "XY": "25", "MT": "26", "M": "26"}
        x = aliases.get(x.upper(), x)
        try:
            return str(int(float(x)))
        except ValueError:
            return x

    bim["chr"] = bim["chr"].map(norm_chr)
    bim["bp"] = pd.to_numeric(bim["bp"]).astype(np.int64)

    # File regioni: 2 o 3 colonne come nello script DiGAS
    regions = []
    with open(region_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.rstrip().split("\t")
            if len(x) == 2:
                x = [x[0] + "-" + x[1], x[0], x[1]]
            if len(x) != 3:
                raise ValueError(f"Riga non valida nel file regioni: {line.rstrip()}")

            chrom1, start = x[1].split(":", 1)
            chrom2, end = x[2].split(":", 1)
            chrom1, chrom2 = norm_chr(chrom1), norm_chr(chrom2)
            if chrom1 != chrom2:
                raise ValueError(f"Regione su cromosomi diversi: {x[0]}")
            regions.append((str(x[0]), chrom1, int(start), int(end)))

    regions = pd.DataFrame(regions, columns=["Gene", "chr", "start", "end"])
    return subjects, bim, regions


def select_snps(sig_file, bim, regions, fold_dir):
    sig = pd.read_csv(sig_file, sep=r"\s+", dtype=str)
    sig["Pvalue"] = pd.to_numeric(sig["Pvalue"], errors="coerce")
    sig = sig[sig["Pvalue"] <= 0.05].copy()

    keep = (
        ((sig["Cat1"] == "CN") & sig["Cat2"].isin(["AD", "MCI"]))
        | ((sig["Cat2"] == "CN") & sig["Cat1"].isin(["AD", "MCI"]))
    )
    sig_cn = sig[keep].copy()
    sig_cn.to_csv(fold_dir / "DiGAS_CN_contrasts.tsv", sep="\t", index=False)

    genes = set(sig_cn["Gene"].dropna().astype(str))
    if not genes:
        raise ValueError("Nessuna regione p<=0.05 nei confronti CN-vs-AD/CN-vs-MCI")

    selected = regions[regions["Gene"].isin(genes)]
    mapping = []
    snps = []
    seen = set()

    for row in selected.itertuples(index=False):
        hit = bim[(bim["chr"] == row.chr) & (bim["bp"] > row.start) & (bim["bp"] < row.end)]
        if hit.empty:
            continue

        rsids = hit["snp"].astype(str).tolist()
        mapping.append({
            "Genes": row.Gene,
            "SNPsCount": len(rsids),
            "SNPsNames": ",".join(f"{row.chr}:{x}" for x in hit["bp"].astype(str)),
            "SNPsNames_Converter": ",".join(rsids),
        })

        for snp in rsids:
            if snp not in seen:
                seen.add(snp)
                snps.append(snp)

    if not snps:
        raise ValueError("Le regioni selezionate non contengono SNP nel BIM")

    pd.DataFrame(mapping).to_csv(
        fold_dir / "selected_regions_SnpToGene.tsv", sep="\t", index=False
    )
    (fold_dir / "selected_snps.txt").write_text("\n".join(snps) + "\n")
    return snps, len(genes), len(sig_cn), len(sig)


def recode_genotypes(root, dataset, plink, snps, pheno_file, work_dir, ad_vs_cn=False):
    work_dir.mkdir(parents=True, exist_ok=True)
    snp_file = work_dir / "snps.txt"
    snp_file.write_text("\n".join(snps) + "\n")
    log = work_dir / "plink.log"

    core = work_dir / "Core"
    raw_prefix = work_dir / "Core_recoded"

    run_command(
        [plink, "--bfile", root / dataset, "--extract", snp_file, "--make-bed", "--out", core],
        log,
    )
    run_command(
        [plink, "--bfile", core, "--recodeA", "--out", raw_prefix],
        log,
    )

    raw = pd.read_csv(str(raw_prefix) + ".raw", sep=r"\s+", dtype={"IID": str})
    pheno = pd.read_csv(
        pheno_file, sep="\t", header=None, names=["IID", "OriginalGroup"], dtype=str
    )
    if ad_vs_cn:
        pheno = pheno[pheno["OriginalGroup"].isin(["CN", "AD"])].copy()
        pheno["Group"] = pheno["OriginalGroup"]
    else:
        pheno["Group"] = pheno["OriginalGroup"].map({"CN": "CN", "MCI": "AD", "AD": "AD"})
    raw = raw.merge(pheno, on="IID", how="inner")

    meta = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "OriginalGroup", "Group"}
    features = [c for c in raw.columns if c not in meta]
    X = raw[features].apply(pd.to_numeric, errors="coerce")
    return X, raw[["IID", "OriginalGroup", "Group"]]


def preprocess(X_train, X_test):
    med = X_train.median(axis=0, skipna=True)
    cols = med[med.notna()].index
    X_train = X_train[cols].fillna(med[cols])
    X_test = X_test[cols].fillna(med[cols])

    bad = X_train.columns[X_train.isna().any() | X_test.isna().any()]
    if len(bad):
        X_train = X_train.drop(columns=bad)
        X_test = X_test.drop(columns=bad)

    if X_train.shape[1] == 0:
        raise ValueError("Nessuna feature utilizzabile dopo l'imputazione")

    return X_train.replace(2, 1), X_test.replace(2, 1)


def tune_knn(X_raw, y, fold):
    n_cv = min(5, int(pd.Series(y).value_counts().min()))
    if n_cv < 2:
        return 1

    cv = StratifiedKFold(n_splits=n_cv, shuffle=True, random_state=fold)
    max_k = min(10, len(y) - math.ceil(len(y) / n_cv))
    scores = []

    for k in range(1, max_k + 1):
        fold_scores = []
        for tr, va in cv.split(X_raw, y):
            Xtr, Xva = preprocess(
                X_raw.iloc[tr].reset_index(drop=True),
                X_raw.iloc[va].reset_index(drop=True),
            )
            clf = KNeighborsClassifier(n_neighbors=k)
            clf.fit(Xtr, y.iloc[tr])
            fold_scores.append(accuracy_score(y.iloc[va], clf.predict(Xva)))
        scores.append(np.mean(fold_scores))

    return int(np.argmax(scores) + 1)


def metrics(y_true, y_pred, prefix):
    acc = accuracy_score(y_true, y_pred)
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=["AD", "CN"], zero_division=0
    )
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=["CN", "AD"]).ravel()

    return {
        f"{prefix}_Accuracy": acc,
        f"{prefix}_BalancedAccuracy": balanced_accuracy_score(y_true, y_pred),
        f"{prefix}_MCC": matthews_corrcoef(y_true, y_pred),
        f"{prefix}_F1": float(np.mean(f1)),
        f"{prefix}_F1-AD": f1[0],
        f"{prefix}_Precision-AD": precision[0],
        f"{prefix}_Recall-AD": recall[0],
        f"{prefix}_F1-CN": f1[1],
        f"{prefix}_Precision-CN": precision[1],
        f"{prefix}_Recall-CN": recall[1],
        f"{prefix}_TP-AD": tp,
        f"{prefix}_TN-AD": tn,
        f"{prefix}_FP-AD": fp,
        f"{prefix}_FN-AD": fn,
    }


def run_fold(
    fold, train_idx, test_idx, subjects, bim, regions, root, dataset, pheno_file,
    region_file, plink, digas_script, digas_iterations, case_dir, total_splits,
    keep_work, mode, ad_vs_cn
):
    label = "fold" if mode == "nested" else "sampling"
    print(f"    {label} {fold}/{total_splits} START")
    fold_dir = case_dir / f"{label}_{fold:02d}"
    fold_dir.mkdir(parents=True, exist_ok=True)

    train = subjects.iloc[train_idx].copy()
    test = subjects.iloc[test_idx].copy()

    # Resume: reuse a fold finished in an earlier (interrupted) run, if it is complete and on the same split.
    pred_name = "oof_predictions.csv" if mode == "nested" else "test_predictions.csv"
    m_file, p_file, t_file = fold_dir / "metrics.csv", fold_dir / pred_name, fold_dir / "outer_test_subjects.tsv"
    if m_file.exists() and p_file.exists() and t_file.exists():
        try:
            old_test = pd.read_csv(t_file, sep="\t", dtype={"IID": str})
            d_old = pd.read_csv(m_file)
            o_old = pd.read_csv(p_file, dtype={"IID": str})
            if (set(old_test.IID) == set(test.IID.astype(str))
                    and len(d_old) == len(METHODS) and (d_old.Status == "OK").all()
                    and len(o_old) == len(METHODS) * len(test)):
                print(f"    {label} {fold}/{total_splits} REUSED (already complete)")
                return fold, d_old.to_dict("records"), o_old.to_dict("records")
        except Exception as exc:
            print(f"    {label} {fold}/{total_splits} cannot reuse ({exc}), recomputing")

    train.to_csv(fold_dir / "outer_train_subjects.tsv", sep="\t", index=False)
    test.to_csv(fold_dir / "outer_test_subjects.tsv", sep="\t", index=False)

    base = {
        "Mode": mode, "Dataset": dataset, "Region": Path(region_file).stem,
        "Target": "AD_vs_CN" if ad_vs_cn else "CN_vs_AD+MCI",
        "Fold": fold, "TotalSplits": total_splits,
        "N_outer_train": len(train), "N_outer_test": len(test),
        "Train_CN_original": int((train.OriginalGroup == "CN").sum()),
        "Train_MCI_original": int((train.OriginalGroup == "MCI").sum()),
        "Train_AD_original": int((train.OriginalGroup == "AD").sum()),
        "Test_CN_original": int((test.OriginalGroup == "CN").sum()),
        "Test_MCI_original": int((test.OriginalGroup == "MCI").sum()),
        "Test_AD_original": int((test.OriginalGroup == "AD").sum()),
    }

    details, oof = [], []
    plink_work = fold_dir / "plink_features"
    # timing: selection measured with GNU time (wall, CPU, peak RSS of the selection process tree);
    # post-selection (SNP mapping + plink recode + preprocessing + 6 classifiers) by wall clock / process CPU
    t_fold, load1 = time.perf_counter(), os.getloadavg()[0]
    sel_time = fold_dir / "selection_time.txt"
    timer = ["/usr/bin/time", "-o", sel_time, "-f", "%e %U %S %M"] if Path("/usr/bin/time").exists() else []

    try:
        # DiGAS soltanto sui soggetti del training fold
        train_pheno = fold_dir / "train_pheno.txt"
        train[["IID", "OriginalGroup"]].to_csv(train_pheno, sep="\t", header=False, index=False)
        out = fold_dir / "digas_train"
        run_command(
            timer + [
                sys.executable, digas_script,
                "--bfile", root / dataset,
                "--region", region_file,
                "--pheno", train_pheno,
                "--out", out,
                "--iteration", str(digas_iterations),
            ],
            fold_dir / "digas_train.log",
            cwd=root,
        )

        sig_file = Path(str(out) + "_significant.txt")
        all_file = Path(str(out) + "_all.txt")
        if not sig_file.exists():
            raise FileNotFoundError(sig_file)

        t_sel_end, cpu_post0 = time.perf_counter(), time.process_time()
        snps, n_regions, n_sig_cn, n_sig_all = select_snps(sig_file, bim, regions, fold_dir)
        X, info = recode_genotypes(root, dataset, plink, snps, pheno_file, plink_work, ad_vs_cn)

        train_ids = set(train.IID.astype(str))
        test_ids = set(test.IID.astype(str))
        train_mask = info.IID.astype(str).isin(train_ids).to_numpy()
        test_mask = info.IID.astype(str).isin(test_ids).to_numpy()

        if train_mask.sum() != len(train_ids) or test_mask.sum() != len(test_ids):
            raise ValueError("Numero di soggetti diverso dopo il merge con i genotipi")

        Xtr_raw = X.loc[train_mask].reset_index(drop=True)
        Xte_raw = X.loc[test_mask].reset_index(drop=True)
        ytr = info.loc[train_mask, "Group"].reset_index(drop=True)
        yte = info.loc[test_mask, "Group"].reset_index(drop=True)
        test_info = info.loc[test_mask].reset_index(drop=True)
        Xtr, Xte = preprocess(Xtr_raw, Xte_raw)

        common = {
            **base,
            "N_significant_rows_all_contrasts": n_sig_all,
            "N_significant_rows_CN_contrasts": n_sig_cn,
            "N_selected_regions": n_regions,
            "N_selected_SNPs": len(snps),
            "N_features_after_preprocessing": Xtr.shape[1],
            "DiGAS_significant_file": str(sig_file),
            "DiGAS_all_file": str(all_file),
        }

        for method in METHODS:
            try:
                k = tune_knn(Xtr_raw, ytr, fold) if method == "KNN" else 5
                if method == "SVM_LINEAR":
                    clf = svm.SVC(kernel="linear")
                elif method == "SVM_POLY":
                    clf = svm.SVC(kernel="poly")
                elif method == "LDA":
                    clf = LDA(n_components=1)
                elif method == "DT":
                    clf = DecisionTreeClassifier(random_state=fold)
                else:
                    clf = KNeighborsClassifier(n_neighbors=k)

                clf.fit(Xtr, ytr)
                pred_tr = clf.predict(Xtr)
                pred_te = clf.predict(Xte)

                row = {
                    **common, "Method": method,
                    "KNN_k": k if method == "KNN" else np.nan,
                    "Status": "OK", "Message": "",
                }
                row.update(metrics(ytr, pred_tr, "Train"))
                row.update(metrics(yte, pred_te, "Test"))
                details.append(row)

                for i, pred in enumerate(pred_te):
                    oof.append({
                        "Mode": mode, "Dataset": dataset, "Region": Path(region_file).stem,
                        "Target": "AD_vs_CN" if ad_vs_cn else "CN_vs_AD+MCI",
                        "Fold": fold, "Method": method,
                        "IID": str(test_info.loc[i, "IID"]),
                        "OriginalGroup": str(test_info.loc[i, "OriginalGroup"]),
                        "TrueLabel": str(yte.iloc[i]), "PredictedLabel": str(pred),
                        "Correct": int(str(yte.iloc[i]) == str(pred)),
                        "KNN_k": k if method == "KNN" else np.nan,
                        "Status": "OK",
                    })
            except Exception as exc:
                details.append({
                    **common, "Method": method, "KNN_k": np.nan,
                    "Status": "ERROR", "Message": str(exc),
                })

        try:
            w = sel_time.read_text().split()[-4:] if sel_time.exists() else ["nan"] * 4
            t_end = time.perf_counter()
            pd.DataFrame([{
                "Fold": fold, "Selection_script": Path(digas_script).name,
                "Selection_wall_s": float(w[0]), "Selection_user_s": float(w[1]), "Selection_sys_s": float(w[2]),
                "Selection_peakRSS_MB": round(float(w[3]) / 1024, 1),
                "PostSelection_wall_s": round(t_end - t_sel_end, 1),
                "PostSelection_python_cpu_s": round(time.process_time() - cpu_post0, 1),
                "Fold_wall_s": round(t_end - t_fold, 1), "Host_load1_at_start": round(load1, 2),
                "N_train": len(train), "N_selected_regions": n_regions, "N_selected_SNPs": len(snps),
            }]).to_csv(fold_dir / "timing.csv", index=False)
        except Exception as exc:
            print(f"    timing not written: {exc}", file=sys.stderr)
        print(f"    {label} {fold}/{total_splits} DONE | regions={n_regions} SNPs={len(snps)}")

    except Exception as exc:
        print(f"    {label} {fold}/{total_splits} ERROR: {exc}", file=sys.stderr)
        for method in METHODS:
            details.append({
                **base, "Method": method, "KNN_k": np.nan,
                "Status": "ERROR", "Message": str(exc),
            })
    finally:
        if details:
            pd.DataFrame(details).to_csv(fold_dir / "metrics.csv", index=False)
        if oof:
            pred_name = "oof_predictions.csv" if mode == "nested" else "test_predictions.csv"
            pd.DataFrame(oof).to_csv(fold_dir / pred_name, index=False)
        if not keep_work:
            shutil.rmtree(plink_work, ignore_errors=True)

    return fold, details, oof


def summarize(details, predictions, total_splits, mode):
    rows = []
    for method in METHODS:
        d = details[(details.Method == method) & (details.Status == "OK")]
        row = {"Method": method}
        if mode == "nested":
            row["OuterFolds_requested"] = total_splits
            row["OuterFolds_successful"] = len(d)
        else:
            row["SamplingRepeats_requested"] = total_splits
            row["SamplingRepeats_successful"] = len(d)

        for col in [
            "Train_Accuracy", "Train_BalancedAccuracy", "Train_MCC",
            "Train_F1", "Train_F1-AD", "Train_F1-CN",
            "Test_Accuracy", "Test_BalancedAccuracy", "Test_MCC",
            "Test_F1", "Test_F1-AD", "Test_F1-CN",
            "N_selected_regions", "N_selected_SNPs", "N_features_after_preprocessing",
        ]:
            values = pd.to_numeric(d[col], errors="coerce").dropna().tolist() if col in d else []
            row[col + "_mean"] = round(statistics.mean(values), 4) if values else np.nan
            row[col + "_stdev"] = round(statistics.stdev(values), 4) if len(values) > 1 else (0.0 if values else np.nan)

        if method == "KNN":
            ks = pd.to_numeric(d.KNN_k, errors="coerce").dropna().astype(int).tolist()
            row["KNN_k_values"] = ",".join(map(str, ks))
            row["KNN_k_median"] = float(np.median(ks)) if ks else np.nan
        else:
            row["KNN_k_values"] = ""
            row["KNN_k_median"] = np.nan

        pp = predictions[(predictions.Method == method) & (predictions.Status == "OK")] if not predictions.empty else pd.DataFrame()
        if not pp.empty:
            prefix = "Pooled_OOF" if mode == "nested" else "Pooled_Test"
            pooled = metrics(pp.TrueLabel.astype(str), pp.PredictedLabel.astype(str), prefix)
            row[prefix + "_N"] = len(pp)
            row[prefix + "_unique_subjects"] = pp.IID.astype(str).nunique()
            row.update({k: round(float(v), 4) for k, v in pooled.items()})
        else:
            prefix = "Pooled_OOF" if mode == "nested" else "Pooled_Test"
            row[prefix + "_N"] = 0
            row[prefix + "_unique_subjects"] = 0

        rows.append(row)
    return pd.DataFrame(rows)

def run_dataset(args, root, dataset, region_file, plink, digas_script, results_dir):
    pheno_file = root / DATASETS[dataset]
    subjects, bim, regions = load_data(
        root, dataset, pheno_file, region_file, ad_vs_cn=args.ad_vs_cn
    )

    if args.sampling_70_30 and args.sampling_90_10:
        raise ValueError("Usare solo uno tra --sampling-70-30 e --sampling-90-10")

    if args.sampling_70_30 or args.sampling_90_10:
        if args.sampling_repeats < 1:
            raise ValueError("--sampling-repeats deve essere >= 1")

        test_size = 0.30 if args.sampling_70_30 else 0.10
        mode = "sampling_70_30" if args.sampling_70_30 else "sampling_90_10"
        work_name = "sampling_work" if args.sampling_70_30 else "sampling90_work"

        idx = np.arange(len(subjects))
        specs = []
        for i in range(args.sampling_repeats):
            tr, te = train_test_split(
                idx, test_size=test_size, stratify=subjects.OriginalGroup,
                random_state=SEED + i,
            )
            specs.append((i + 1, (tr, te)))
        total_splits = args.sampling_repeats
    else:
        mode = "nested"
        counts = subjects.OriginalGroup.value_counts()
        if counts.min() < args.outer_folds:
            raise ValueError(f"Impossibile fare {args.outer_folds} fold: {counts.to_dict()}")
        cv = StratifiedKFold(n_splits=args.outer_folds, shuffle=True, random_state=SEED)
        specs = list(enumerate(
            cv.split(np.zeros(len(subjects)), subjects.OriginalGroup), start=1
        ))
        total_splits = args.outer_folds
        work_name = "nested_work"

    target_tag = "_ADvsCN" if args.ad_vs_cn else ""
    case_dir = results_dir / work_name / f"{dataset}_{region_file.stem}{target_tag}"
    case_dir.mkdir(parents=True, exist_ok=True)

    jobs = min(args.outer_jobs, total_splits)
    if jobs < 1:
        raise ValueError("--outer-jobs deve essere >= 1")

    results = []
    if jobs == 1:
        for fold, (tr, te) in specs:
            results.append(run_fold(
                fold, tr, te, subjects, bim, regions, root, dataset, pheno_file,
                region_file, plink, digas_script, args.digas_iterations,
                case_dir, total_splits, args.keep_work, mode, args.ad_vs_cn,
            ))
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = [
                ex.submit(
                    run_fold, fold, tr, te, subjects, bim, regions, root, dataset,
                    pheno_file, region_file, plink, digas_script, args.digas_iterations,
                    case_dir, total_splits, args.keep_work, mode, args.ad_vs_cn,
                )
                for fold, (tr, te) in specs
            ]
            for f in as_completed(futures):
                results.append(f.result())

    detail_rows, prediction_rows = [], []
    for _, d, o in sorted(results, key=lambda x: x[0]):
        detail_rows.extend(d)
        prediction_rows.extend(o)

    details = pd.DataFrame(detail_rows)
    predictions = pd.DataFrame(prediction_rows)
    if not details.empty:
        details = details.sort_values(["Fold", "Method"]).reset_index(drop=True)
    if not predictions.empty:
        predictions = predictions.sort_values(["Method", "Fold", "IID"]).reset_index(drop=True)

    return summarize(details, predictions, total_splits, mode), details, predictions, mode

def main():
    args = parse_args()
    root = Path(args.root).expanduser().resolve()

    region_file = Path(args.region_file)
    if not region_file.is_absolute():
        region_file = (root / region_file).resolve()

    plink = Path(args.plink)
    if not plink.is_absolute():
        plink = (root / str(plink).removeprefix("./")).resolve()

    digas_script = Path(args.digas_script)
    if not digas_script.is_absolute():
        digas_script = (root / str(digas_script).removeprefix("./")).resolve()

    results_dir = Path(args.results_dir)
    if not results_dir.is_absolute():
        results_dir = (root / results_dir).resolve()

    required = [region_file, plink, digas_script]
    for dataset in args.datasets:
        required += [root / f"{dataset}.bed", root / f"{dataset}.bim", root / f"{dataset}.fam", root / DATASETS[dataset]]
    missing = [str(x) for x in required if not Path(x).exists()]
    if missing:
        raise FileNotFoundError("File mancanti:\n" + "\n".join(missing))

    if args.sampling_70_30 and args.sampling_90_10:
        raise ValueError("Usare solo uno tra --sampling-70-30 e --sampling-90-10")
    if args.sampling_70_30:
        mode_dir = "sampling_70_30"
    elif args.sampling_90_10:
        mode_dir = "sampling_90_10"
    else:
        mode_dir = "nested"
    (results_dir / mode_dir).mkdir(parents=True, exist_ok=True)
    (results_dir / "details" / mode_dir).mkdir(parents=True, exist_ok=True)

    all_summary = []
    manifest = []

    for dataset in args.datasets:
        if args.sampling_70_30:
            mode = "sampling_70_30"
        elif args.sampling_90_10:
            mode = "sampling_90_10"
        else:
            mode = "nested"
        target = "AD_vs_CN" if args.ad_vs_cn else "CN_vs_AD+MCI"
        print(f"\n[{mode.upper()}] {dataset} | {region_file.stem} | {target}")
        row = {
            "Mode": mode, "Target": target, "Dataset": dataset, "Region": region_file.stem,
            "RegionFile": str(region_file), "OuterFolds": args.outer_folds,
            "SamplingRepeats": args.sampling_repeats,
            "RandomState": SEED, "DiGAS_iterations": args.digas_iterations,
            "OuterJobs": args.outer_jobs, "Status": "STARTED", "Message": "",
        }

        try:
            summary, details, predictions, mode = run_dataset(
                args, root, dataset, region_file, plink, digas_script, results_dir
            )
            target_tag = "_ADvsCN" if args.ad_vs_cn else ""
            tag = f"{dataset}_{region_file.stem}{target_tag}"
            summary.to_csv(results_dir / mode / f"{tag}.csv", index=False)
            detail_suffix = "folds" if mode == "nested" else "splits"
            details.to_csv(
                results_dir / "details" / mode / f"{tag}_{detail_suffix}.csv", index=False
            )
            pred_suffix = "oof_predictions" if mode == "nested" else "test_predictions"
            predictions.to_csv(
                results_dir / "details" / mode / f"{tag}_{pred_suffix}.csv", index=False
            )

            x = summary.copy()
            x.insert(0, "Region", region_file.stem)
            x.insert(0, "Dataset", dataset)
            all_summary.append(x)
            row["Status"] = "OK"
        except Exception as exc:
            row["Status"] = "ERROR"
            row["Message"] = str(exc)
            print(f"ERROR: {exc}", file=sys.stderr)

        manifest.append(row)
        pd.DataFrame(manifest).to_csv(results_dir / "run_manifest.csv", index=False)

    if all_summary:
        if args.sampling_70_30:
            out_name = "ALL_SAMPLING_70_30.csv"
        elif args.sampling_90_10:
            out_name = "ALL_SAMPLING_90_10.csv"
        else:
            out_name = "ALL_NESTED.csv"
        pd.concat(all_summary, ignore_index=True).to_csv(results_dir / out_name, index=False)

    print(f"\nRisultati: {results_dir}")


if __name__ == "__main__":
    main()
