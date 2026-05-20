# -*- coding: utf-8 -*-
"""
East Asian case-control PRS evaluation using logistic regression.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2


BASE_DIR = Path("/path/to/eas_prs_evaluation")
case_control_path = BASE_DIR / "eas_case_control_long.csv"
prs_path = BASE_DIR / "EAS_5e-8_col.sscore"
output_path = BASE_DIR / "PRS_R2_EAS_5e-8.csv"

cancer_prs_lst = [
    "bladder",
    "breast",
    "cervical",
    "colorectal",
    "endometrial",
    "kidney",
    "lung",
    "melanoma",
    "oral",
    "ovarian",
    "prostate",
    "thyroid",
    "meta",
    "pan",
]
cancer_outcome_lst = [
    "All cancer",
    "Liver cancer",
    "Lung cancer",
    "Breast cancer",
    "Kidney cancer",
    "Thyroid cancer",
    "Gastric cancer",
    "Pancreatic cancer",
    "Colorectal cancer",
    "Oral cancer",
]
female_cancers = ["Breast cancer", "Ovarian cancer", "Endometrial cancer", "Cervical cancer"]
male_cancers = ["Prostate cancer"]


def nagelkerke_delta_r2(ll0, llc, llp, n):
    denom = 1 - np.exp((2 / n) * ll0)
    cs_full = 1 - np.exp((2 / n) * (ll0 - llp))
    cs_cov = 1 - np.exp((2 / n) * (ll0 - llc))
    nk_full = cs_full / denom
    nk_cov = cs_cov / denom
    return nk_full - nk_cov


phenotype = pd.read_csv(case_control_path)
prs_score = pd.read_csv(prs_path, sep="\t")

score_cols = dict(zip(cancer_prs_lst, [f"SCORE{i}_AVG" for i in range(1, len(cancer_prs_lst) + 1)]))
prs_keep = ["IID"]
for prs_name in cancer_prs_lst:
    named_avg_col = f"{prs_name}_AVG"
    score_avg_col = score_cols[prs_name]
    avg_col = named_avg_col if named_avg_col in prs_score.columns else score_avg_col
    prs_score[prs_name] = prs_score["ALLELE_CT"] * prs_score[avg_col]
    prs_keep.append(prs_name)

result = []
for prs_name in cancer_prs_lst:
    for outcome in cancer_outcome_lst:
        df = phenotype.loc[phenotype["Cancer"] == outcome].copy()
        df = pd.merge(df, prs_score[["IID", prs_name]], on="IID", how="inner")

        if outcome in female_cancers:
            df = df[df["SEX"] == 1].copy()
            covars = ["age", "PC1", "PC2", "PC3", "PC4", "PC5"]
        elif outcome in male_cancers:
            df = df[df["SEX"] == 0].copy()
            covars = ["age", "PC1", "PC2", "PC3", "PC4", "PC5"]
        else:
            covars = ["age", "SEX", "PC1", "PC2", "PC3", "PC4", "PC5"]

        df = df.dropna(subset=["case", prs_name] + covars).copy()
        df["case"] = df["case"].astype(int)
        df["PRS_sd"] = (df[prs_name] - df[prs_name].mean()) / df[prs_name].std(ddof=0)
        df["constant"] = 1

        n = len(df)
        n_case = int(df["case"].sum())
        if n == 0 or n_case == 0 or n_case == n:
            result.append([prs_name, outcome, n, n_case, np.nan, np.nan, np.nan])
            continue

        try:
            null_model = sm.Logit(df["case"], np.ones(len(df))).fit(disp=0)
            model_covariate = sm.Logit(df["case"], df[["constant"] + covars]).fit(disp=0)
            model_prs = sm.Logit(df["case"], df[["constant", "PRS_sd"] + covars]).fit(disp=0)

            ll0 = null_model.llf
            llc = model_covariate.llf
            llp = model_prs.llf

            lr_stat = 2 * (llp - llc)
            p_lr = float(chi2.sf(lr_stat, df=1))
            delta_nk = nagelkerke_delta_r2(ll0, llc, llp, n)

            mcfadden_full = 1 - (llp / ll0)
            mcfadden_cov = 1 - (llc / ll0)
            delta_mcfadden = mcfadden_full - mcfadden_cov

            result.append([prs_name, outcome, n, n_case, delta_nk * 100, delta_mcfadden * 100, p_lr])
        except Exception:
            result.append([prs_name, outcome, n, n_case, np.nan, np.nan, np.nan])

result_df = pd.DataFrame(
    result,
    columns=["PRS", "Cancer", "N", "Events", "Nagelkerke_R2_pct", "McFadden_R2_pct", "P_LR"],
)
result_df.to_csv(output_path, index=False)
