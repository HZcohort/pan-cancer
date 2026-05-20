from pathlib import Path

import pandas as pd


UKB_RS = Path("/path/to/ukb_rs.csv")
PRS_DIR = Path("/path/to/PRS_EraSOR")
OUT_FILE = PRS_DIR / "score_merged.txt.gz"

ukb_rs = pd.read_csv(UKB_RS)
rs_ukb_dict = dict(ukb_rs[["rs", "name"]].values)

score_merged = []
for cancer in [
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
]:
    try:
        score = pd.read_csv(PRS_DIR / cancer / f"{cancer}_sbrc.txt", sep="\t")
    except FileNotFoundError:
        print(cancer)
        continue
    score["SNP"] = score["SNP"].map(rs_ukb_dict)
    score = score[["SNP", "A1", "BETA"]].dropna()
    score.rename(columns={"BETA": cancer}, inplace=True)
    try:
        score_merged = pd.merge(score_merged, score, on=["SNP", "A1"], how="outer")
    except TypeError:
        score_merged = score

score_merged.to_csv(OUT_FILE, sep="\t", index=False)
