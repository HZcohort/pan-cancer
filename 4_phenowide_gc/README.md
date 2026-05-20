# Phenome-Wide Genetic Correlation

This folder contains the LDSC workflow used to estimate genetic correlations between the latent pan-cancer factor and publicly available disease or quantitative-trait GWAS summary statistics.

## Scripts

- `phewas_gc.sh` munges the pan-cancer GWAS summary statistics, runs LDSC genetic-correlation analysis against all trait summary statistics in the configured phenome-wide GWAS directory, and calls the merge script.
- `ldsc_merge.py` parses LDSC log output and combines genetic-correlation estimates, standard errors, p-values, heritability estimates, lambda values, and intercepts into one table.

## Expected Inputs

- Pan-cancer GWAS summary statistics with columns `snp`, `p`, `A1`, `A2`, `N`, and `beta`.
- HapMap3 SNP list for LDSC munging.
- LDSC software directory containing `munge_sumstats.py` and `ldsc.py`.
- European LDSC LD-score and regression-weight files.
- A directory of munged phenome-wide GWAS summary statistics named as `*.sumstats.gz`.

