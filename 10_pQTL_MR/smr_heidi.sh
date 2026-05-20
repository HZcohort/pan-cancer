#!/usr/bin/env bash
set -euo pipefail

SMR="/path/to/smr"
REF_BFILE="/path/to/1000G_EUR/g1000_EUR"
OUT_DIR="/path/to/pqtl_mr/smr_heidi"
THREADS="8"
INPUT_TABLE="/path/to/smr_input_manifest.tsv"

mkdir -p "${OUT_DIR}"

while IFS=$'\t' read -r dataset trait gwas_ma besd_prefix
do
    if [[ "${dataset}" == "dataset" ]]; then
        continue
    fi

    "${SMR}" \
        --bfile "${REF_BFILE}" \
        --gwas-summary "${gwas_ma}" \
        --beqtl-summary "${besd_prefix}" \
        --out "${OUT_DIR}/SMR_${dataset}_${trait}" \
        --thread-num "${THREADS}"
done < "${INPUT_TABLE}"
