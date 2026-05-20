# Risk Loci Annotation

This folder contains the scripts used for post-GWAS processing of the pan-cancer factor GWAS. These scripts correspond to the manuscript section on identification of LD-independent risk loci and GWAS Catalog novelty assessment.

## Script Order

1. `LiftRsNumber.py`
   - Standardizes the pan-cancer GWAS summary statistics.
   - Updates merged rs identifiers using dbSNP history files.
   - Aligns variants to a 1000 Genomes European reference panel.

2. External PLINK clumping
   - Identifies genome-wide significant index variants and LD clumps from the processed pan-cancer GWAS.
   - Uses the manuscript clumping parameters: index variants at `p-value < 5e-8`, secondary variants at `p-value <= 1e-5`, LD threshold `r2 >= 0.1`, and a 2 Mb window.
   - Produces `pan_locus.clumped` and `pan_indsig.clumped` files for downstream steps.

3. `loci_identify.py`
   - Converts PLINK clumps into genomic intervals.
   - Merges nearby or overlapping clump intervals.

4. `loci_verify.py`
   - Selects the lead variant with the smallest p-value within each locus.
   - Produces verified lead-variant tables, BED files, and `loci.parameter` files.

5. External LocusZoom plotting
   - Uses the verified `loci.parameter` file to generate regional association plots.

6. `catalog.py`
   - Queries GWAS Catalog associations near each pan-cancer locus.
   - Uses hg19-to-hg38 liftover before matching against GWAS Catalog coordinates.

## Expected Inputs

- Pan-cancer GWAS summary statistics with SNP ID, alleles, effect size, standard error, p-value, and sample-size columns.
- dbSNP rs history files preprocessed as `rs_history.npy` and `rs_merge.npy`.
- 1000 Genomes European PLINK reference files, including a `.bim` file.
- PLINK clumping outputs named with the same prefix used for the processed GWAS.
- GWAS Catalog association file in tab-delimited format with `P-VALUE`, `CHR_ID`, and `CHR_POS` columns.
- Liftover chain files available to `pyliftover`.
