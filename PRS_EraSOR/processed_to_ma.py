from pathlib import Path

import pandas as pd


PROJECT_DIR = Path("/path/to/PRS_EraSOR")
EUR_AFREQ = Path("/path/to/EUR.afreq.gz")

sample_size_adjusted_dict = {
    "melanoma": 190241,
    "cervical": 46708,
    "endometrial": 40250,
    "oral": 48240,
    "ovarian": 104336,
    "thyroid": 32070,
    "prostate": 383422,
    "colorectal": 341185,
    "lung": 192344,
    "bladder": 66001,
    "kidney": 119798,
    "breast": 430306,
    "meta": 1994902,
}

cancer_lst = [
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
]

af = pd.read_csv(EUR_AFREQ, sep="\t")
flip_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}

for cancer in cancer_lst:
    df = pd.read_csv(PROJECT_DIR / f"{cancer}_adjusted.txt", sep="\t")
    if cancer == "meta":
        df = df.rename(columns={"MarkerName": "snp",
                                "Allele1": "A1",
                                "Allele2": "A2",
                                "Effect": "beta",
                                "StdErr": "se",
                                "P-value": "p"})
        df["A1"] = df["A1"].str.upper()
        df["A2"] = df["A2"].str.upper()
    df = pd.merge(df, af[["ID", "REF", "ALT", "ALT_FREQS"]], left_on="snp", right_on="ID", how="left")
    cond1 = df[(df["A1"] == df["ALT"]) | (df["A2"] == df["REF"])].index
    cond2 = df[(df["A1"] == df["REF"]) | (df["A2"] == df["ALT"])].index
    cond3 = df[(df["A1"] == df["ALT"].map(flip_dict)) | (df["A2"] == df["REF"].map(flip_dict))].index
    cond4 = df[(df["A1"] == df["REF"].map(flip_dict)) | (df["A2"] == df["ALT"].map(flip_dict))].index
    df.loc[cond1, "freq"] = df.loc[cond1, "ALT_FREQS"]
    df.loc[cond2, "freq"] = 1 - df.loc[cond2, "ALT_FREQS"]
    df.loc[cond3, "freq"] = df.loc[cond3, "ALT_FREQS"]
    df.loc[cond4, "freq"] = 1 - df.loc[cond4, "ALT_FREQS"]
    df = df.rename(columns={"snp": "SNP", "beta": "b"})
    df["N"] = sample_size_adjusted_dict[cancer]
    df[["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]].to_csv(
        PROJECT_DIR / f"{cancer}.ma", sep="\t", index=False
    )
