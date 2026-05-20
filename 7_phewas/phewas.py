# %% PheWAS of the pan-cancer PRS
from pathlib import Path

import numpy as np
import pandas as pd
import DiNetxify as dnt


BASE_DIR = Path("/path/to/phewas_analysis")
OUTPUT_DIR = BASE_DIR / "result"
DINETXIFY_DIR = BASE_DIR / "dinetxify_data"

BASELINE_PATH = BASE_DIR / "baseline.csv"
PRS_PATH = BASE_DIR / "pan_eras_PRS.sscore"
UNRELATED_IDS_PATH = BASE_DIR / "unrelated_european_ids.csv"
HOSPITAL_RECORDS_PATH = BASE_DIR / "hospital_inpatient_death_registry_icd10.csv"
CANCER_RECORDS_PATH = BASE_DIR / "cancer_registry_icd10.csv"

PHECODE_VERSION = "1.3a"
N_PROCESS = 25
N_CASES_EXPOSED_THRESHOLD = 100

COVARIATES = [
    "age",
    "social_cate",
    "BMI",
    "income",
    "smoking",
    "drinking",
    "physical",
    "diet",
    "pca1",
    "pca2",
    "pca3",
    "pca4",
    "pca5",
]

SYSTEMS_EXCLUDED = [
    "congenital anomalies",
    "symptoms",
    "others",
    "pregnancy complications",
    "injuries & poisonings",
]


def read_unrelated_ids(path):
    unrelated = pd.read_csv(path)
    for col in ["IID", "eid", "participant_id"]:
        if col in unrelated.columns:
            return set(unrelated[col].astype(str))
    return set(unrelated.iloc[:, 0].astype(str))


def prepare_phenotype():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DINETXIFY_DIR.mkdir(parents=True, exist_ok=True)

    baseline = pd.read_csv(BASELINE_PATH)
    prs = pd.read_csv(PRS_PATH, sep="\t")

    if {"BETA_AVG", "ALLELE_CT"}.issubset(prs.columns):
        prs["PRS"] = prs["BETA_AVG"] * prs["ALLELE_CT"]
    elif "PRS" not in prs.columns:
        score_cols = [col for col in prs.columns if col.endswith("_AVG")]
        if len(score_cols) != 1:
            raise ValueError("PRS file must contain PRS, BETA_AVG, or one *_AVG score column.")
        prs["PRS"] = prs[score_cols[0]] * prs["ALLELE_CT"]

    prs = prs.dropna(subset=["PRS"]).copy()
    q3 = prs["PRS"].quantile(0.75)
    prs["case"] = np.where(prs["PRS"] >= q3, 1, 0)

    unrelated_ids = read_unrelated_ids(UNRELATED_IDS_PATH)
    prs["IID"] = prs["IID"].astype(str)
    prs_remain = prs[prs["IID"].isin(unrelated_ids)].copy()

    baseline["eid"] = baseline["eid"].astype(str)
    phenotype = pd.merge(
        prs_remain[["IID", "case"]],
        baseline,
        left_on="IID",
        right_on="eid",
        how="inner",
    )

    phenotype_columns = [
        "eid",
        "case",
        "sex",
        "age",
        "social_cate",
        "BMI",
        "income",
        "smoking",
        "drinking",
        "physical",
        "diet",
        "pca1",
        "pca2",
        "pca3",
        "pca4",
        "pca5",
        "date_attend",
        "date_end",
    ]
    phenotype[phenotype_columns].to_csv(OUTPUT_DIR / "phenotype_pan.csv", index=False)


def run_phewas(cancer="pan"):
    col_dict = {
        "Participant ID": "eid",
        "Exposure": "case",
        "Sex": "sex",
        "Index date": "date_attend",
        "End date": "date_end",
    }

    data = dnt.DiseaseNetworkData(
        study_design="cohort",
        phecode_level=1,
        date_fmt="%Y-%m-%d",
        phecode_version=PHECODE_VERSION,
    )
    data.phenotype_data(
        phenotype_data_path=str(OUTPUT_DIR / f"phenotype_{cancer}.csv"),
        column_names=col_dict,
        covariates=COVARIATES,
    )

    data.merge_medical_records(
        medical_records_data_path=str(HOSPITAL_RECORDS_PATH),
        diagnosis_code="ICD-10-WHO",
        column_names={
            "Participant ID": "eid",
            "Diagnosis code": "diag_icd10",
            "Date of diagnosis": "dia_date",
        },
    )
    data.merge_medical_records(
        medical_records_data_path=str(CANCER_RECORDS_PATH),
        diagnosis_code="ICD-10-WHO",
        column_names={
            "Participant ID": "eid",
            "Diagnosis code": "ICD10",
            "Date of diagnosis": "date",
        },
    )
    data.save(str(DINETXIFY_DIR / cancer))

    data = dnt.DiseaseNetworkData()
    data.load(str(DINETXIFY_DIR / cancer))
    data.modify_phecode_level(2)

    phewas_result = dnt.phewas(
        data=data,
        covariates=[
            "age",
            "pca1",
            "pca2",
            "pca3",
            "pca4",
            "pca5",
            "sex",
            "social_cate",
            "BMI",
            "income",
            "smoking",
            "drinking",
            "physical",
            "diet",
        ],
        n_threshold=N_CASES_EXPOSED_THRESHOLD,
        n_process=N_PROCESS,
        correction="fdr_bh",
        cutoff=0.05,
        system_exl=SYSTEMS_EXCLUDED,
        log_file=str(OUTPUT_DIR / "phewas.log"),
        lifelines_disable=True,
    )
    phewas_result.to_csv(OUTPUT_DIR / f"phewas_all_{cancer}.csv", index=False)

    phewas_eligible = phewas_result[
        phewas_result["N_cases_exposed"] >= N_CASES_EXPOSED_THRESHOLD
    ].copy()
    phewas_eligible.to_csv(OUTPUT_DIR / f"phewas_eligible_{cancer}.csv", index=False)


if __name__ == "__main__":
    prepare_phenotype()
    run_phewas("pan")
