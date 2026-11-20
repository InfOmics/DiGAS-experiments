#!/usr/bin/env bash
# DiGAS genotype QC with PLINK 1.9.
#
# Usage:
#   bash scripts/qc/qc_plink19.sh <input_bfile_prefix> <output_bfile_prefix>
# Example (from repository root):
#   bash scripts/qc/qc_plink19.sh ADNI3 ADNI3_QC
#
# Requires <input_bfile_prefix>.bed/.bim/.fam.
# Applies all four filters in one PLINK invocation.
# Does NOT perform ancestry filtering or genome-build conversion.
# Note: this is a reproducible QC implementation based on the documented
#       parameters. The exact historical commands for all cohorts have not
#       been independently verified.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input_bfile_prefix> <output_bfile_prefix>" >&2
  exit 2
fi

INPUT="$1"
OUTPUT="$2"
PLINK_BIN="${PLINK_BIN:-./plink}"

for ext in bed bim fam; do
  if [[ ! -f "${INPUT}.${ext}" ]]; then
    echo "Error: missing input file: ${INPUT}.${ext}" >&2
    exit 1
  fi
done

if [[ "$PLINK_BIN" == */* ]]; then
  if [[ ! -x "$PLINK_BIN" ]]; then
    echo "Error: PLINK executable not found or not executable: $PLINK_BIN" >&2
    exit 1
  fi
elif ! command -v "$PLINK_BIN" >/dev/null 2>&1; then
  echo "Error: PLINK executable not found in PATH: $PLINK_BIN" >&2
  exit 1
fi

if [[ "$INPUT" == "$OUTPUT" ]]; then
  echo "Error: input and output prefixes must be different." >&2
  exit 1
fi

for ext in bed bim fam; do
  if [[ -e "${OUTPUT}.${ext}" ]]; then
    echo "Error: output already exists (${OUTPUT}.${ext}); choose another prefix to avoid overwriting data." >&2
    exit 1
  fi
done

mkdir -p "$(dirname "$OUTPUT")"

echo "Input PLINK prefix: $INPUT"
echo "Output PLINK prefix: $OUTPUT"
echo "Applying QC in one PLINK 1.9 run: --geno 0.2 --mind 0.1 --maf 0.05 --hwe 1e-6"

"$PLINK_BIN" \
  --bfile "$INPUT" \
  --geno 0.2 \
  --mind 0.1 \
  --maf 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --out "$OUTPUT"

echo "QC complete. Review ${OUTPUT}.log and the output .bed/.bim/.fam files."
