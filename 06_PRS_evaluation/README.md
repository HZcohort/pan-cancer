# PRS Evaluation

This folder contains the scripts used to evaluate polygenic risk score (PRS) performance in the UK Biobank and in the independent East Asian case-control validation dataset.

## Workflow

`baseline_process.py` prepares the UK Biobank baseline covariate table used in the prospective PRS analysis. It derives age at recruitment, lifestyle covariates, genotype principal components, and censoring dates, and writes `baseline.csv`.

`PRS.py` first constructs wide cancer outcome files from linked hospital inpatient, death registry, and cancer registry ICD-10 records. It then evaluates the pan-cancer, cross-cancer, and site-specific PRSs using Cox proportional hazards models. For each PRS-outcome pair, it reports hazard ratios per standard deviation increase in PRS and the incremental Royston-Sauerbrei R2 from nested model comparison.

`PRS_eas.py` evaluates PRSs in the East Asian validation dataset using logistic regression. It assumes that cohort-specific diagnosis processing has already been completed and uses a prepared long-format case-control table for regression analysis.

## Expected Inputs

The UK Biobank scripts expect:

- a UK Biobank baseline table containing `eid` and raw baseline field columns used in `baseline_process.py`;
- a participant inclusion file with currently valid UK Biobank participant IDs;
- lifestyle variables containing `eid`, `alcohol_gday`, and `smoking_num`;
- genotype principal components containing `eid` and the first five PCs;
- a death or censoring-date table;
- hospital inpatient (death registry data incorporated) and cancer registry ICD-10 records with participant IDs, ICD-10 diagnosis codes, and diagnosis dates;
- `cancer_ICD.xlsx`, containing cancer names and ICD-10 definitions;
- PLINK2 `.sscore` files for the site-specific, cross-cancer, and pan-cancer PRSs;
- an unrelated European-ancestry participant ID list for the final analytic sample.

The East Asian validation script expects:

- `eas_case_control_long.csv`, with columns `IID`, `Cancer`, `case`, `SEX`, `age`, and `PC1` to `PC5`;
- an EAS PLINK2 `.sscore` file containing `IID`, `ALLELE_CT`, and either named PRS columns such as `pan_AVG` or score columns such as `SCORE14_AVG`.

