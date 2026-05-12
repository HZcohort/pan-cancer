# %% data generation
from pathlib import Path
import pandas as pd


BASE_DIR = Path("/path/to/ukb_prs_evaluation")
hospital_inpatient_path = BASE_DIR / "input" / "hospital_inpatient_death_registry_icd10.csv"
cancer_registry_path = BASE_DIR / "input" / "cancer_registry_icd10.csv"
baseline_path = BASE_DIR / "baseline.csv"
icd_definition_path = str(BASE_DIR / "cancer_ICD.xlsx")
output_dir = str(BASE_DIR / "icd10_case_control")

baseline_id_col = "eid"

hospital_id_col = "eid"
hospital_code_col = "diag_icd10"
hospital_date_col = "dia_date"

registry_id_col = "eid"
registry_code_col = "ICD10"
registry_date_col = "date"


def clean_icd10(x):
    return str(x).upper().replace(".", "").strip()


def get_icd_sets(icd_rule):
    icd3 = set()
    icd4 = set()

    for item in str(icd_rule).replace("–", "-").replace("—", "-").split(","):
        item = item.strip()
        if item == "":
            continue

        if "-" in item:
            start, end = item.split("-", 1)
            start = clean_icd10(start)
            end = clean_icd10(end)

            if len(start) == 3 and len(end) == 3:
                for i in range(int(start[1:]), int(end[1:]) + 1):
                    icd3.add(start[0] + str(i).zfill(2))

            elif len(start) == 4 and len(end) == 4:
                for i in range(int(start[3]), int(end[3]) + 1):
                    icd4.add(start[:3] + str(i))

        else:
            code = clean_icd10(item)
            if len(code) == 3:
                icd3.add(code)
            elif len(code) >= 4:
                icd4.add(code[:4])

    return icd3, icd4


if __name__ == "__main__":
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cancer_icd = pd.read_excel(icd_definition_path)
    cancer_icd = cancer_icd[["Cancer type", "ICD-10 codes"]].copy()

    baseline = pd.read_csv(baseline_path)
    baseline = baseline.rename(columns={baseline_id_col: "eid"})

    hospital = pd.read_csv(hospital_inpatient_path)
    hospital = hospital[[hospital_id_col, hospital_code_col, hospital_date_col]].copy()
    hospital.columns = ["eid", "icd10", "diag_date"]

    registry = pd.read_csv(cancer_registry_path)
    registry = registry[[registry_id_col, registry_code_col, registry_date_col]].copy()
    registry.columns = ["eid", "icd10", "diag_date"]

    records = pd.concat([hospital, registry], ignore_index=True)
    records["icd10"] = records["icd10"].map(clean_icd10)
    records["icd10_3"] = records["icd10"].str[:3]
    records["icd10_4"] = records["icd10"].str[:4]
    records["diag_date"] = pd.to_datetime(records["diag_date"])

    cancer_wide = baseline[["eid"]].copy()
    cancer_date_wide = baseline[["eid"]].copy()

    for _, row in cancer_icd.iterrows():
        cancer_name = row["Cancer type"]
        icd_rule = row["ICD-10 codes"]
        cancer_icd3, cancer_icd4 = get_icd_sets(icd_rule)

        case_records = records.loc[
            records["icd10_3"].isin(cancer_icd3) | records["icd10_4"].isin(cancer_icd4),
            ["eid", "diag_date"],
        ].copy()
        case_dates = (
            case_records.groupby("eid", as_index=False)["diag_date"]
            .min()
            .rename(columns={"diag_date": f"{cancer_name}_date"})
        )

        cancer_date_wide = cancer_date_wide.merge(case_dates, on="eid", how="left")
        cancer_wide[cancer_name] = cancer_date_wide[f"{cancer_name}_date"].notna().astype(int)

        print(cancer_name, cancer_wide[cancer_name].sum())

    for col in cancer_date_wide.columns:
        if col != "eid":
            cancer_wide[col] = cancer_date_wide[col]

    cancer_wide.to_csv(output_dir / "cancer_case_control_wide.csv", index=False)

    print("Output folder:", output_dir)
    print("Totally cancer-free controls:", (cancer_wide["All cancer"] == 0).sum())

# %% PRS -- prospective Cox models with incremental Royston-Sauerbrei R2
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from statsmodels.api import PHReg


BASE_DIR = Path("/path/to/ukb_prs_evaluation")
cancer_path = BASE_DIR / "icd10_case_control" / "cancer_case_control_wide.csv"
baseline_path = BASE_DIR / "baseline.csv"
prs_path = BASE_DIR / "PRS_erase.sscore"
pan_prs_path = BASE_DIR / "pan_eras_PRS.sscore"
unrelated_ids_path = BASE_DIR / "unrelated_european_ids.csv"
output_path = BASE_DIR / "PRS_R_2.csv"

baseline_cols = ["eid", "sex", "date_attend", "date_end", "age", "pca1", "pca2", "pca3", "pca4", "pca5"]
female_cancers = ["Breast cancer", "Endometrial cancer", "Ovarian cancer", "Cervical cancer"]
male_cancers = ["Prostate cancer"]
prs_avg_cols = {
    "meta": "meta_AVG",
    "bladder": "bladder_AVG",
    "cervical": "cervical_AVG",
    "colorectal": "colorectal_AVG",
    "breast": "breast_AVG",
    "endometrial": "endometrial_AVG",
    "lung": "lung_AVG",
    "kidney": "kidney_AVG",
    "melanoma": "melanoma_AVG",
    "ovarian": "ovarian_AVG",
    "prostate": "prostate_AVG",
    "thyroid": "thyroid_AVG",
    "oral": "oral_AVG",
}


def royston_r2_d(time, event, linear_predictor):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    linear_predictor = np.asarray(linear_predictor, dtype=float)
    n = len(linear_predictor)
    z = norm.ppf((np.arange(1, n + 1) - 3 / 8) / (n + 0.25))
    order = np.argsort(linear_predictor, kind="mergesort")
    qhat = np.empty(n, dtype=float)
    sorted_lp = linear_predictor[order]

    start = 0
    while start < n:
        end = start + 1
        while end < n and sorted_lp[end] == sorted_lp[start]:
            end += 1
        qhat[order[start:end]] = z[start:end].mean()
        start = end

    royston_df = pd.DataFrame({"time": time, "event": event, "qhat": qhat})
    royston_fit = PHReg.from_formula("time ~ qhat", status="event", data=royston_df, ties="breslow").fit(disp=0)
    d_stat = float(royston_fit.params[0] * np.sqrt(8 / np.pi))
    r2_d = float(d_stat**2 / (np.pi**2 / 6 + d_stat**2))
    return d_stat, r2_d


df_cancer = pd.read_csv(cancer_path)
cancer_cols = [col for col in df_cancer.columns if col != "eid" and not col.endswith("_date")]
baseline = pd.read_csv(baseline_path)
df = pd.merge(baseline[baseline_cols], df_cancer, on="eid", how="inner")
df["date_attend"] = pd.to_datetime(df["date_attend"])
df["date_end"] = pd.to_datetime(df["date_end"])
for cancer in cancer_cols:
    df[f"{cancer}_date"] = pd.to_datetime(df[f"{cancer}_date"])

# read site-specific PRS
prs = pd.read_csv(prs_path, sep="\t")
prs["eid"] = prs["IID"]
prs_keep = ["eid"]
for prs_name, avg_col in prs_avg_cols.items():
    prs_col = f"PRS_{prs_name}"
    prs[prs_col] = prs[avg_col] * prs["ALLELE_CT"]
    prs[prs_col] = (prs[prs_col] - prs[prs_col].mean()) / prs[prs_col].std(ddof=0)
    prs_keep.append(prs_col)
df = pd.merge(df, prs[prs_keep], on="eid", how="inner")

# read pan PRS
prs_pan = pd.read_csv(pan_prs_path, sep="\t")
prs_pan["eid"] = prs_pan["IID"]
prs_pan["PRS_pan"] = prs_pan["BETA_AVG"] * prs_pan["ALLELE_CT"]
prs_pan["PRS_pan"] = (prs_pan["PRS_pan"] - prs_pan["PRS_pan"].mean()) / prs_pan["PRS_pan"].std(ddof=0)
df = pd.merge(df, prs_pan[["eid", "PRS_pan"]], on="eid", how="inner")

# keep the unrelated European-ancestry analysis sample
id_included = pd.read_csv(unrelated_ids_path)
df = df[df["eid"].isin(id_included["eid"])].copy()

prs_names = ["pan"] + list(prs_avg_cols.keys())
results = []

for cancer_outcome in cancer_cols:
    print(cancer_outcome)

    date_col = f"{cancer_outcome}_date"
    prs_model = df.copy()

    prevalent = prs_model[date_col].notna() & (prs_model[date_col] <= prs_model["date_attend"])
    prs_model = prs_model.loc[~prevalent].copy()

    prs_model["event"] = (
        prs_model[date_col].notna()
        & (prs_model[date_col] > prs_model["date_attend"])
        & (prs_model[date_col] <= prs_model["date_end"])
    ).astype(int)
    prs_model["end_date"] = prs_model["date_end"]
    prs_model.loc[prs_model["event"] == 1, "end_date"] = prs_model.loc[prs_model["event"] == 1, date_col]
    prs_model["time"] = ((prs_model["end_date"] - prs_model["date_attend"]).dt.days + 1) / 365.25

    covariates = ["age", "sex", "pca1", "pca2", "pca3", "pca4", "pca5"]
    if cancer_outcome in male_cancers:
        prs_model = prs_model[prs_model["sex"] == 0].copy()
        covariates = ["age", "pca1", "pca2", "pca3", "pca4", "pca5"]
    if cancer_outcome in female_cancers:
        prs_model = prs_model[prs_model["sex"] == 1].copy()
        covariates = ["age", "pca1", "pca2", "pca3", "pca4", "pca5"]

    prs_model = prs_model.dropna(subset=["time", "event"] + covariates)
    prs_model = prs_model[prs_model["time"] > 0].copy()

    for prs_name in prs_names:
        prs_col = f"PRS_{prs_name}"
        model_df = prs_model.dropna(subset=[prs_col]).copy()
        n = len(model_df)
        n_event = int(model_df["event"].sum())

        if n == 0 or n_event == 0 or n_event == n:
            results.append(
                [prs_name, cancer_outcome, n, n_event, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
            )
            continue

        try:
            cov_formula = "time ~ " + " + ".join(covariates)
            full_formula = "time ~ " + " + ".join([prs_col] + covariates)

            model_cov = PHReg.from_formula(cov_formula, status="event", data=model_df, ties="breslow").fit(disp=0)
            model_full = PHReg.from_formula(full_formula, status="event", data=model_df, ties="breslow").fit(disp=0)

            lp_cov = model_cov.model.exog.dot(model_cov.params)
            lp_full = model_full.model.exog.dot(model_full.params)

            _, r2_cov = royston_r2_d(model_df["time"], model_df["event"], lp_cov)
            _, r2_full = royston_r2_d(model_df["time"], model_df["event"], lp_full)

            lr_stat = 2 * (model_full.llf - model_cov.llf)
            p_lr = float(chi2.sf(lr_stat, df=1))

            prs_idx = model_full.model.exog_names.index(prs_col)
            beta = float(model_full.params[prs_idx])
            se = float(model_full.bse[prs_idx])
            hr = float(np.exp(beta))
            hr_l95 = float(np.exp(beta - 1.96 * se))
            hr_u95 = float(np.exp(beta + 1.96 * se))
            p_prs = float(model_full.pvalues[prs_idx])

            results.append(
                [
                    prs_name,
                    cancer_outcome,
                    n,
                    n_event,
                    r2_cov * 100,
                    r2_full * 100,
                    (r2_full - r2_cov) * 100,
                    hr,
                    hr_l95,
                    hr_u95,
                    p_prs,
                    p_lr,
                ]
            )

        except Exception:
            results.append(
                [prs_name, cancer_outcome, n, n_event, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
            )


df_result = pd.DataFrame(
    results,
    columns=[
        "PRS",
        "Cancer",
        "N",
        "Events",
        "Royston_Sauerbrei_R2_covariates_pct",
        "Royston_Sauerbrei_R2_full_pct",
        "Delta_Royston_Sauerbrei_R2_pct",
        "HR_per_SD",
        "HR_L95",
        "HR_U95",
        "P_PRS",
        "P_LR",
    ],
)
df_result.to_csv(output_path, index=False)

