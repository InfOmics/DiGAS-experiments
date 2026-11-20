#!/usr/bin/env python3
"""Summarise the DiGAS vs SKAT comparison (results/tool_comparison/<arm>_<tool>_<region>/).

Primary measures (the ones reported in the paper, used for the comparison): accuracy, micro-F1, macro-F1.
Computed from the out-of-fold predictions: pooled over all subjects, plus mean/sd and the 10 per-fold values.
Micro-F1 is computed independently from the confusion counts and checked against accuracy.
Secondary measures (shown, NOT used for the comparison): balanced accuracy, MCC.
Reference: majority-class accuracy of the same subjects.

  .venv/bin/python scripts/analysis/compare_tools.py      -> results/tool_comparison/summary_*.csv + printout
"""
import re
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score, matthews_corrcoef

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tool_comparison"
JOB = re.compile(r"^(leaky|nested)_(DiGAS|SKAT)_(.+)$")


def micro_f1_counts(t, p):
    """Micro-F1 from summed per-class TP/FP/FN (independent of accuracy_score)."""
    tp = fp = fn = 0
    for c in np.union1d(t, p):
        tp += np.sum((p == c) & (t == c)); fp += np.sum((p == c) & (t != c)); fn += np.sum((p != c) & (t == c))
    return 2 * tp / (2 * tp + fp + fn)


def scores(t, p):
    t, p = np.asarray(t), np.asarray(p)
    acc = accuracy_score(t, p); mic = micro_f1_counts(t, p)
    assert abs(acc - mic) < 1e-12, f"micro-F1 {mic} != accuracy {acc}"
    return dict(acc=acc, microF1=mic, macroF1=f1_score(t, p, average="macro"),
                bal_acc=balanced_accuracy_score(t, p), MCC=matthews_corrcoef(t, p))


def main():
    rows, folds = [], []
    for d in sorted(OUT.iterdir()):
        m = JOB.match(d.name)
        if not (d.is_dir() and m and (OUT / f"{d.name}.done").exists()):
            continue
        arm, tool, region = m.groups()
        for f in (d / "details" / "nested").glob("*_oof_predictions.csv"):
            P = pd.read_csv(f)
            P = P[P.Status == "OK"]
            for meth, g in P.groupby("Method"):
                maj = g.TrueLabel.value_counts(normalize=True).max()
                s = scores(g.TrueLabel, g.PredictedLabel)
                fs = pd.DataFrame([dict(fold=k, **scores(h.TrueLabel, h.PredictedLabel))
                                   for k, h in g.groupby("Fold")])
                for _, r in fs.iterrows():
                    folds.append(dict(region=region, arm=arm, tool=tool, method=meth, **r))
                rows.append(dict(region=region, arm=arm, tool=tool, method=meth, n=len(g), majority_acc=maj,
                                 **{f"pooled_{k}": v for k, v in s.items()},
                                 **{f"fold_{k}_mean": fs[k].mean() for k in s},
                                 **{f"fold_{k}_sd": fs[k].std(ddof=1) for k in s},
                                 fold_acc_values=" ".join(f"{x:.3f}" for x in fs.acc),
                                 fold_macroF1_values=" ".join(f"{x:.3f}" for x in fs.macroF1)))
    if not rows:
        print("no finished jobs yet"); return
    S = pd.DataFrame(rows).sort_values(["region", "method", "arm", "tool"])
    F = pd.DataFrame(folds)
    S.to_csv(OUT / "summary_all.csv", index=False); F.to_csv(OUT / "summary_folds.csv", index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
    print("PRIMARY (comparison): pooled out-of-fold accuracy / micro-F1 / macro-F1; fold mean±sd; majority-class acc")
    P = S[["region", "method", "arm", "tool", "majority_acc", "pooled_acc", "pooled_microF1", "pooled_macroF1",
           "fold_acc_mean", "fold_acc_sd", "fold_macroF1_mean", "fold_macroF1_sd"]]
    print(P.round(3).to_string(index=False))
    print("\nSECONDARY (not used for the comparison): pooled balanced accuracy / MCC")
    print(S[["region", "method", "arm", "tool", "pooled_bal_acc", "pooled_MCC"]].round(3).to_string(index=False))
    # paired clean DiGAS vs clean SKAT, per region x method, on the primary measures (fold-paired: same folds)
    out = []
    for (reg, meth), g in F[F.arm == "nested"].groupby(["region", "method"]):
        a = g[g.tool == "DiGAS"].set_index("fold"); b = g[g.tool == "SKAT"].set_index("fold")
        common = a.index.intersection(b.index)
        if len(common) < 2:
            continue
        for k in ("acc", "macroF1"):
            dlt = a.loc[common, k] - b.loc[common, k]
            from scipy.stats import wilcoxon
            p = wilcoxon(dlt).pvalue if np.any(dlt != 0) else 1.0
            out.append(dict(region=reg, method=meth, measure=k, DiGAS_minus_SKAT_mean=dlt.mean(),
                            folds_DiGAS_better=int((dlt > 0).sum()), folds_SKAT_better=int((dlt < 0).sum()),
                            ties=int((dlt == 0).sum()), wilcoxon_p=p))
    if out:
        C = pd.DataFrame(out); C.to_csv(OUT / "summary_paired_nested.csv", index=False)
        print("\nCLEAN DiGAS vs clean SKAT, fold-paired (accuracy, macro-F1)")
        print(C.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
