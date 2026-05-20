# Pathway Enrichment

This folder contains the public workflow scripts for pathway enrichment analyses of the pan-cancer factor GWAS and the prioritized pan-cancer risk genes.

## Scripts

- `magma_gene_set_analysis.sh`: runs MAGMA SNP-to-gene analysis for the pan-cancer GWAS, performs MAGMA gene-set enrichment against MSigDB pathway annotations, and calls the merge script.
- `MAGMA_geneset_merge.py`: merges MAGMA gene-set outputs from canonical pathways, chemical and genetic perturbation signatures, hallmark pathways, and Gene Ontology collections; adds MSigDB metadata; and calculates collection-specific q-values.
- `gene_based_enrichment.py`: performs the second-stage mapped-gene enrichment analysis using GSEApy, restricted to pathways significant in the MAGMA gene-set screen.

## Expected Inputs

The MAGMA workflow expects a pan-cancer GWAS input file with `snp` and `p` columns, a European LD reference panel, MAGMA gene annotation files, MAGMA-formatted MSigDB gene-set annotations, and MSigDB JSON metadata files.

The mapped-gene enrichment script expects the merged MAGMA gene-set table (`gene_set_MAGMA.csv`), a consensus mapped-gene table with a `Gene` column, an Ensembl-to-gene-symbol table with a `NAME` column used as the background gene universe, and MSigDB JSON files for GO, canonical pathways, and hallmark gene sets.

