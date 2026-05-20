# pQTL-MR

This folder contains the analysis scripts for the proteome-wide Mendelian randomization analysis of blood proteins against the pan-cancer factor GWAS.

The public workflow starts from harmonized cis-pQTL/GWAS input files. Each harmonized input file should be a gzipped tab-delimited table with one row per variant-protein association and the following columns: `feature_unique`, `chr`, `bp`, `A1_qtl`, `A2_qtl`, `A1F`, `beta_qtl`, `se_qtl`, `p_qtl`, `snp`, `A1`, `A2`, `beta`, `se`, and `p`. Optional columns such as `N_gwas`, `N_xQTL`, and `type_gwas` can be retained.

Scripts:

- `initial_filter.R` performs LD clumping within each protein using `ieugwasr::ld_clump`, the specified pQTL p-value threshold, and a 1000 Genomes European reference panel.
- `prepare_mr_ld_inputs.sh` connects the clumping step with GCTA genotype extraction, producing the filtered `.gz`, SNP allele `.tsv`, and `.xmat.gz` files used by the MR script.
- `ivw_mr_parallel.R` performs generalized least-squares inverse-variance weighted MR using the `MendelianRandomization` package while accounting for LD correlation among retained instruments.
- `select_replicated_proteins.py` applies FDR correction to IVW results and identifies proteins reaching significance in at least two pQTL datasets with consistent effect direction.
- `smr_heidi.sh` runs SMR/HEIDI verification for candidate proteins from prepared SMR `.ma` and BESD inputs.

Expected reference inputs include PLINK/GCTA-compatible 1000 Genomes European reference files, a population allele-frequency reference used during harmonization, harmonized cis-pQTL/GWAS tables for AGES, deCODE, Fenland, and UK Biobank pQTL datasets, and prepared SMR `.ma` and BESD files for SMR/HEIDI verification.
