# PheWAS of the Pan-Cancer PRS

This folder contains the script used to conduct the phenome-wide association study (PheWAS) of the EraSOR-adjusted pan-cancer polygenic risk score (PRS) in the UK Biobank analytic cohort.

## Workflow

`phewas.py` prepares the exposure file by merging the pan-cancer PRS with the same baseline covariate table and unrelated European-ancestry participant list used in the PRS evaluation. Participants in the top quartile of the pan-cancer PRS distribution are coded as the high genetic risk group, and the remaining participants are used as the reference group.

The script then uses `DiNetxify` to map linked ICD-10 healthcare records to phecodes, construct outcome-specific incident disease cohorts, and run PheWAS models adjusted for demographic, lifestyle, and genetic principal-component covariates. Phecodes from congenital anomalies, symptoms, pregnancy complications, injuries and poisonings, and other non-specific categories are excluded. The final eligible results are restricted to phecodes with at least 100 incident cases in the high genetic risk group.

## Expected Inputs

The script expects:

- `baseline.csv`, containing `eid`, `sex`, `age`, `social_cate`, `BMI`, `income`, `smoking`, `drinking`, `physical`, `diet`, `pca1` to `pca5`, `date_attend`, and `date_end`;
- `pan_eras_PRS.sscore`, a PLINK2 score file containing `IID`, `ALLELE_CT`, and either `BETA_AVG`, `PRS`, or one `*_AVG` PRS column;
- `unrelated_european_ids.csv`, containing participant IDs for the unrelated European-ancestry analytic cohort;
- hospital inpatient and death-registry ICD-10 records with participant ID, diagnosis code, and diagnosis date columns;
- cancer-registry ICD-10 records with participant ID, diagnosis code, and diagnosis date columns.

