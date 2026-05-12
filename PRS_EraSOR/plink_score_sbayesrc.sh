#!/bin/bash
#SBATCH -p cn
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -J PRS_score
#SBATCH -o logs/PRS_score.o
#SBATCH -e logs/PRS_score.e

set -euo pipefail

PLINK2="/path/to/plink2"
UKB_PGEN_PREFIX="/path/to/ukb_imp_qc"
UKB_FREQ="/path/to/ukb_imp_qc.acount"
SCORE_FILE="/path/to/PRS_EraSOR/score_merged.txt.gz"
OUT_PREFIX="/path/to/output/PRS_erase"

"${PLINK2}" \
  --bpfile "${UKB_PGEN_PREFIX}" \
  --out "${OUT_PREFIX}" \
  --read-freq "${UKB_FREQ}" \
  --score "${SCORE_FILE}" 1 2 no-mean-imputation header-read \
  --score-col-nums 3-16 \
  --threads 40
