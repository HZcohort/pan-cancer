import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


WORK_DIR = Path("/path/to/erasor_workdir")
BASE_GWAS_DIR = Path("/path/to/base_gwas")
UKB_GWAS_DIR = Path("/path/to/ukb_gwas")
ADJUSTED_DIR = Path("/path/to/erasor_adjusted_gwas")
ERASOR_SCRIPT = Path("/path/to/EraSOR.py")
EUR_AFREQ = Path("/path/to/EUR.afreq.gz")
LDSC_REF_DIR = Path("/path/to/ld_score/LDSC/EUR")
HM3_SNPLIST = Path("/path/to/w_hm3.snplist")
PAN_FACTOR_GWAS = Path("/path/to/pan_eras.gz")


UKB_CASE_COUNTS = {
    "melanoma": 14255,
    "cervical": 1032,
    "endometrial": 5908,
    "oral": 624,
    "ovarian": 3808,
    "thyroid": 1516,
    "prostate": 23165,
    "colorectal": 12883,
    "lung": 8480,
    "bladder": 3048,
    "kidney": 3952,
}

ADJUSTED_SAMPLE_SIZES = {
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
}

ADJUSTED_CANCERS = [
    "melanoma",
    "cervical",
    "endometrial",
    "oral",
    "ovarian",
    "thyroid",
    "prostate",
    "colorectal",
    "lung",
    "bladder",
    "kidney",
]


def add_sample_size_to_ukb_targets():
    for cancer in ADJUSTED_CANCERS:
        path = UKB_GWAS_DIR / f"{cancer}_UKB.gz"
        df = pd.read_csv(path, sep="\t")
        df["N"] = UKB_CASE_COUNTS[cancer]
        df.to_csv(path, index=False, sep="\t")


def run_erasor_adjustment():
    for cancer in ADJUSTED_CANCERS:
        cmd = [
            "python",
            str(ERASOR_SCRIPT),
            "--base", str(BASE_GWAS_DIR / f"{cancer}.gz"),
            "--base-snp", "snp",
            "--base-signed-sumstats", "beta,0",
            "--base-a1", "A1",
            "--base-a2", "A2",
            "--base-p", "p",
            "--base-N-col", "N",
            "--target", str(UKB_GWAS_DIR / f"{cancer}_UKB.gz"),
            "--target-snp", "variant_id",
            "--target-signed-sumstats", "beta,0",
            "--target-a1", "effect_allele",
            "--target-a2", "other_allele",
            "--target-p", "p_value",
            "--target-N-col", "N",
            "--ref-ld-chr", str(LDSC_REF_DIR) + "/",
            "--w-ld-chr", str(LDSC_REF_DIR) + "/",
            "--out", str(ADJUSTED_DIR / f"{cancer}_adjusted"),
        ]
        subprocess.run(cmd, check=True)


def add_adjusted_sample_size():
    for cancer in ADJUSTED_SAMPLE_SIZES:
        path = ADJUSTED_DIR / f"{cancer}_adjusted.gz"
        df = pd.read_csv(path, sep="\t")
        df["N"] = ADJUSTED_SAMPLE_SIZES[cancer]
        df.to_csv(path, index=False, sep="\t")


def generate_beta_se_files():
    freq = pd.read_csv(EUR_AFREQ, sep="\t")
    freq_dict = dict(freq[["ID", "ALT_FREQS"]].values)

    for cancer in ADJUSTED_CANCERS:
        print(cancer)
        df = pd.read_csv(ADJUSTED_DIR / f"{cancer}_adjusted.gz", sep="\t")
        df["maf"] = df["SNP"].map(freq_dict)
        df["deno"] = 2 * df["maf"] * (1 - df["maf"]) * (df["N"] + df["Z"] ** 2)
        df["deno"] = np.sqrt(df["deno"])
        df["beta"] = df["Z"] / df["deno"]
        df["se"] = df["beta"] / df["Z"]
        df = df.rename(columns={"SNP": "snp", "P": "p"})
        df[["snp", "beta", "se", "A1", "A2", "p"]].to_csv(
            ADJUSTED_DIR / f"{cancer}_adjusted.txt", index=False, sep="\t"
        )


def create_pan_factor_ma():
    af = pd.read_csv(EUR_AFREQ, sep="\t")
    df = pd.read_csv(PAN_FACTOR_GWAS, sep="\t")
    df = pd.merge(df, af[["ID", "ALT", "ALT_FREQS"]], left_on="snp", right_on="ID", how="left")
    index_a1_alt = df[df["A1"] == df["ALT"]].index
    df.loc[index_a1_alt, "freq"] = df.loc[index_a1_alt, "ALT_FREQS"]
    df.loc[~df.index.isin(index_a1_alt), "freq"] = 1 - df.loc[~df.index.isin(index_a1_alt), "ALT_FREQS"]
    df = df.rename(columns={"snp": "SNP", "beta": "b"})
    df["N"] = 6127
    df[["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]].to_csv(
        WORK_DIR / "pan.ma", sep="\t", index=False
    )


def main():
    add_sample_size_to_ukb_targets()
    run_erasor_adjustment()
    add_adjusted_sample_size()
    generate_beta_se_files()
    create_pan_factor_ma()


if __name__ == "__main__":
    main()
