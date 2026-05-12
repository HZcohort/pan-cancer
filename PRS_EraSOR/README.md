# PRS EraSOR Workflow

This folder contains the scripts used to construct EraSOR-adjusted polygenic risk scores for the revised PRS analyses.

## Workflow

- `eras_process.py` adds UK Biobank target sample sizes, applies EraSOR to site-specific GWAS summary statistics that included UK Biobank participants, converts EraSOR Z-score output back to beta/SE format, and prepares the pan-cancer factor GWAS in SBayesRC `.ma` format.
- `gSEM.R` reconstructs the pan-cancer factor GWAS from the EraSOR-adjusted cancer-specific summary statistics.
- `metal_meta_analysis.metal` and `run_metal.sh` perform standard-error-mode METAL meta-analysis across the 12 PRS-ready cancer GWAS to generate the cross-cancer GWAS used for the comparison PRS.
- `processed_to_ma.py` converts the 12 site-specific GWAS and the cross-cancer meta-analysis GWAS into SBayesRC `.ma` format.
- `SBayesRC.sh` estimates PRS weights for the 12 site-specific PRSs, the cross-cancer PRS, and the pan-cancer PRS.
- `sbayesrc_to_plink_score.py` merges SBayesRC output weights into a multi-score PLINK input file.
- `plink_score_sbayesrc.sh` calculates individual-level SBayesRC PRSs.
- `pt_pipeline_eas.py` constructs clumping-and-thresholding PRS score files for East Asian validation at the genome-wide significance threshold.
- `convert.py` maps rsID-based P+T score files to UK Biobank variant IDs when PLINK scoring requires UK Biobank variant names.

## Expected Inputs

- Base GWAS summary statistics for the 12 constituent cancers.
- UK Biobank cancer GWAS summary statistics for cancers requiring EraSOR adjustment.
- Breast cancer GWAS summary statistics formatted consistently with the EraSOR-adjusted outputs.
- Pan-cancer factor GWAS reconstructed from EraSOR-adjusted cancer-specific summary statistics.
- European allele-frequency reference file and LDSC European LD-score resources.
- East Asian 1000 Genomes PLINK reference panel for P+T validation PRS construction.
- SBayesRC LD reference files and baseline annotation file.
- UK Biobank rsID-to-variant-name mapping file and UK Biobank imputed genotype files for PLINK scoring.
