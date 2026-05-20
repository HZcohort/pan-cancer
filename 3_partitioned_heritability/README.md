# Partitioned Heritability

This folder contains the stratified LD score regression workflow used to estimate SNP-heritability enrichment for the latent pan-cancer factor.

## Script

- `partitioned_ldsc.sh` munges the pan-cancer GWAS summary statistics and runs LDSC enrichment analyses for the baselineLD model, tissue and cell-type-specific gene-expression annotations, and cell-type-specific histone-mark annotations.

## Expected Inputs

- Pan-cancer GWAS summary statistics with columns `snp`, `p`, `A1`, `A2`, `N`, and `beta`.
- HapMap3 SNP list for LDSC munging.
- LDSC software directory containing `munge_sumstats.py` and `ldsc.py`.
- European LD-score resources for baselineLD, baseline, regression weights, multi-tissue expression annotations, and multi-tissue chromatin annotations.

