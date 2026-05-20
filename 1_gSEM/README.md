# GenomicSEM Pan-Cancer Factor Analysis

This folder contains the scripts used for the gSEM analysis described in the manuscript sections **Deriving the latent pan-cancer factor using genomic structural equation modelling** and **GWAS of Pan-cancer latent factor from genomic structural equation modelling**.

## Script Order

1. `main_gsem.R`
   - Estimates the LDSC genetic covariance and sampling covariance matrices.
   - Processes cancer GWAS summary statistics for GenomicSEM.
   - Runs chromosome-specific multivariate GWAS for the selected one-factor model.

2. `model_eval.R`
   - Fits the LDSC-based one-, two-, and three-factor models.
   - Prints parameter estimates and model-fit statistics.

3. `hdl_compare.R`
   - Estimates the HDL covariance structure.
   - Fits the HDL-based one-factor model used as a sensitivity comparison.

## Expected Inputs

- GenomicSEM munged summary statistics named:
  `bladder.sumstats.gz`, `breast.sumstats.gz`, `cervical.sumstats.gz`, `colorectal.sumstats.gz`, `endometrial.sumstats.gz`, `kidney.sumstats.gz`, `lung.sumstats.gz`, `melanoma.sumstats.gz`, `oral.sumstats.gz`, `ovarian.sumstats.gz`, `prostate.sumstats.gz`, and `thyroid.sumstats.gz`.
- Processed GWAS files named:
  `bladder.txt`, `breast.txt`, `cervical.txt`, `colorectal.txt`, `endometrial.txt`, `kidney.txt`, `lung.txt`, `melanoma.txt`, `oral.txt`, `ovarian.txt`, `prostate.txt`, and `thyroid.txt`.
- European LDSC LD-score and weight files.
- HDL LD reference files for the sensitivity analysis.
- A 1000 Genomes European reference file for `GenomicSEM::sumstats()`.

