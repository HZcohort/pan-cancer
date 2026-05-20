#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUT_DIR="/path/to/MAGMA_enrich"
REF_BFILE="/path/to/reference/g1000_EUR"
PAN_GWAS_PVAL="/path/to/pan_cancer_magma_input"
GENE_ANNOT="/path/to/MAGMA/EUR.annot.genes.annot"
MAGMA_GENE_SET_DIR="/path/to/MAGMA/gene_set_annotations"
MSIGDB_JSON_DIR="/path/to/MSigDB/json"
N_CASES="6127"

mkdir -p "${OUT_DIR}/gtex" "${OUT_DIR}/gene_set" "${OUT_DIR}/data"

magma \
  --bfile "${REF_BFILE}" \
  --pval "${PAN_GWAS_PVAL}" use=snp,p N="${N_CASES}" \
  --gene-annot "${GENE_ANNOT}" \
  --gene-model snp-wise=mean \
  --out "${OUT_DIR}/data/pan"

magma \
  --gene-results "${OUT_DIR}/data/pan.genes.raw" \
  --set-annot "${MAGMA_GENE_SET_DIR}/cp.out" \
  --out "${OUT_DIR}/gene_set/pan_cp"

magma \
  --gene-results "${OUT_DIR}/data/pan.genes.raw" \
  --set-annot "${MAGMA_GENE_SET_DIR}/cgp.out" \
  --out "${OUT_DIR}/gene_set/pan_cgp"

magma \
  --gene-results "${OUT_DIR}/data/pan.genes.raw" \
  --set-annot "${MAGMA_GENE_SET_DIR}/hallmark.out" \
  --out "${OUT_DIR}/gene_set/pan_hallmark"

magma \
  --gene-results "${OUT_DIR}/data/pan.genes.raw" \
  --set-annot "${MAGMA_GENE_SET_DIR}/GO.out" \
  --out "${OUT_DIR}/gene_set/pan_GO"

python "${SCRIPT_DIR}/MAGMA_geneset_merge.py" \
  --input_path "${OUT_DIR}/gene_set" \
  --output_path "${OUT_DIR}" \
  --trait_list pan \
  --geneset_path "${MSIGDB_JSON_DIR}"
