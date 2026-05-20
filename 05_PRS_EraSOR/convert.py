from pathlib import Path

import pandas as pd


UKB_RS = Path("/path/to/ukb_rs.csv")
SCORE_DIR = Path("/path/to/pt_scores")
OUT_DIR = Path("/path/to/plink_score_files")

ukb_rs = pd.read_csv(UKB_RS)
rs_ukb_dict = dict(ukb_rs[["rs", "name"]].values)

score_files = ["EUR_5e-8.score"]
for file in score_files:
    score = pd.read_csv(SCORE_DIR / file, sep="\t")
    score["SNP"] = score["SNP"].map(rs_ukb_dict)
    score = score.dropna(subset=["SNP"])
    score.to_csv(OUT_DIR / file, sep="\t", index=False)
