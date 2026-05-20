import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_dataset_file(values):
    mapping = {}
    for value in values:
        if "=" not in value:
            raise ValueError("Dataset inputs must use the format DATASET=/path/to/file.csv")
        dataset, path = value.split("=", 1)
        mapping[dataset] = Path(path)
    return mapping


def bh_fdr(pvalues):
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    p_valid = p.loc[valid]
    order = np.argsort(p_valid.values)
    ranked = p_valid.values[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out.loc[p_valid.index[order]] = adjusted
    return out


def standardize_protein_name(value):
    return re.sub(r"_\d+$", "", str(value).strip())


def read_ivw_table(dataset, path, trait):
    df = pd.read_csv(path)
    protein_col = next(
        (col for col in ["Protein", "Feature", "Gene", "feature_unique"] if col in df.columns),
        None,
    )
    if protein_col is None:
        raise ValueError(f"{path} does not contain a protein identifier column")
    required = {"Trait", "IVW_Beta", "IVW_SE", "IVW_P", "IVW_N_Total"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {','.join(sorted(missing))}")
    df = df[df["Trait"].astype(str) == trait].copy()
    df["Dataset"] = dataset
    df["Protein"] = df[protein_col].map(standardize_protein_name)
    df["IVW_P"] = pd.to_numeric(df["IVW_P"], errors="coerce")
    df["IVW_Beta"] = pd.to_numeric(df["IVW_Beta"], errors="coerce")
    df["IVW_SE"] = pd.to_numeric(df["IVW_SE"], errors="coerce")
    df["IVW_N_Total"] = pd.to_numeric(df["IVW_N_Total"], errors="coerce")
    df = df.sort_values("IVW_P", na_position="last").drop_duplicates(["Dataset", "Protein"])
    df["IVW_FDR"] = bh_fdr(df["IVW_P"])
    return df[["Dataset", "Protein", "Trait", "IVW_Beta", "IVW_SE", "IVW_P", "IVW_FDR", "IVW_N_Total"]]


def summarize(dataset_level, fdr_threshold, min_datasets):
    dataset_level["Significant"] = dataset_level["IVW_FDR"] < fdr_threshold
    dataset_level["Direction"] = np.sign(dataset_level["IVW_Beta"])
    summary = (
        dataset_level.groupby("Protein", as_index=False)
        .agg(
            N_datasets_tested=("IVW_P", lambda x: x.notna().sum()),
            N_datasets_significant=("Significant", "sum"),
            Directions=("Direction", lambda x: ",".join(map(str, sorted(set(x.dropna().astype(int)))))),
        )
    )
    direction_sets = dataset_level.dropna(subset=["Direction"]).groupby("Protein")["Direction"].apply(
        lambda x: set(x.astype(int))
    )
    summary["Direction_consistent"] = summary["Protein"].map(
        lambda protein: len(direction_sets.get(protein, set())) <= 1
    )
    summary["Replicated_IVW_candidate"] = (
        (summary["N_datasets_significant"] >= min_datasets) & summary["Direction_consistent"]
    )
    return summary.sort_values(
        ["Replicated_IVW_candidate", "N_datasets_significant", "Protein"],
        ascending=[False, False, True],
    )


def main():
    parser = argparse.ArgumentParser(description="Select replicated pQTL-MR protein associations")
    parser.add_argument("--ivw_results", nargs="+", required=True, help="DATASET=/path/to/IVW_dataset.csv")
    parser.add_argument("--trait", default="pan")
    parser.add_argument("--fdr_threshold", type=float, default=0.05)
    parser.add_argument("--min_datasets", type=int, default=2)
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    dataset_files = parse_dataset_file(args.ivw_results)
    dataset_level = pd.concat(
        [read_ivw_table(dataset, path, args.trait) for dataset, path in dataset_files.items()],
        ignore_index=True,
    )
    summary = summarize(dataset_level, args.fdr_threshold, args.min_datasets)

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    dataset_level.to_csv(output_prefix.with_name(output_prefix.name + "_dataset_level.csv"), index=False)
    summary.to_csv(output_prefix.with_name(output_prefix.name + "_protein_summary.csv"), index=False)


if __name__ == "__main__":
    main()
