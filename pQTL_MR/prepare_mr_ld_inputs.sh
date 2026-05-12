#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="/path/to/harmonized_pqtl_gwas"
OUTPUT_DIR="/path/to/pqtl_mr/data"
PLINK="/path/to/plink"
GCTA="/path/to/gcta64"
REF_BFILE="/path/to/1000G_EUR/g1000_EUR"
R2_THRESH="0.05"
PVAL_THRESH="5e-8"
N_CORES="20"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${OUTPUT_DIR}"

while IFS=$'\t' read -r dataset trait part input_file
do
    if [[ "${dataset}" == "dataset" ]]; then
        continue
    fi

    prefix="${OUTPUT_DIR}/${dataset}.${trait}.Part${part}"

    Rscript "${SCRIPT_DIR}/initial_filter.R" \
        "${input_file}" \
        "${prefix}.gz" \
        "feature_unique" \
        "${R2_THRESH}" \
        "${PVAL_THRESH}" \
        "${PLINK}" \
        "${REF_BFILE}"

    zcat "${prefix}.gz" | awk 'NR>1{print $10,$11}' > "${prefix}.tsv"

    "${GCTA}" \
        --bfile "${REF_BFILE}" \
        --extract "${prefix}.tsv" \
        --update-ref-allele "${prefix}.tsv" \
        --recode \
        --thread-num "${N_CORES}" \
        --out "${prefix}"
done < "/path/to/harmonized_input_manifest.tsv"
