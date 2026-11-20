#!/usr/bin/env python3
"""Leaky protocol through the unchanged nested pipeline: a stand-in for the DiGAS script that ignores the fold's
training subjects and returns a region selection computed beforehand on ALL subjects.

  PRECOMPUTED_PREFIX=<prefix of <prefix>_significant.txt / <prefix>_all.txt> \
  python scripts/nested/digas_classification_nested_10fold.py --digas-script scripts/analysis/precomputed_selection.py ...

Every outer fold then uses the same whole-cohort selection (test subjects included): the protocol of the paper.
Same command line as scripts/nested/digas3-linux_validation.py; only --out is used.
"""
import argparse, os, shutil

p = argparse.ArgumentParser()
for a in ("--bfile", "--region", "--pheno", "--out"):
    p.add_argument(a)
p.add_argument("--iteration"); p.add_argument("--saveTemporary", action="store_false")
args = p.parse_args()
src = os.environ["PRECOMPUTED_PREFIX"]
for suf in ("_significant.txt", "_all.txt"):
    if os.path.exists(src + suf):
        shutil.copyfile(src + suf, args.out + suf)
