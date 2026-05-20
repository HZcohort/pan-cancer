#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
REF_PREFIX = Path("/path/to/Reference_panel/EAS/g1000_EAS")
REF_BIM = REF_PREFIX.with_suffix(".bim")
PLINK2 = Path("/path/to/plink2")
WORK_DIR = BASE_DIR / "pt_work_EAS"
CLUMP_DIR = BASE_DIR / "pt_clump_EAS"
SCORE_DIR = BASE_DIR / "pt_scores_EAS"

MAX_P = 5e-8
THRESHOLDS = [("5e-8", 5e-8)]
TRAITS = [
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


def parse_float(value: str) -> float | None:
    value = value.strip()
    if value == "" or value.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def open_gwas(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def row_get(row: dict[str, str], aliases: tuple[str, ...]) -> str:
    for alias in aliases:
        if alias in row:
            return row[alias]
    return ""


def iter_gwas(trait: str):
    path = BASE_DIR / ("meta.ma" if trait == "meta" else f"{trait}.gz")
    with open_gwas(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        for row in reader:
            row = {key.strip(): value.strip() for key, value in row.items() if key is not None}
            snp = row_get(row, ("SNP", "snp"))
            a1 = row_get(row, ("A1", "a1")).upper()
            a2 = row_get(row, ("A2", "a2")).upper()
            beta = parse_float(row_get(row, ("b", "beta", "BETA")))
            p_value = parse_float(row_get(row, ("p", "P", "Pval_Estimate")))
            if p_value is None or beta is None or p_value > MAX_P:
                continue
            if not snp or len(a1) != 1 or len(a2) != 1:
                continue
            yield snp, a1, a2, beta, p_value


def ref_rsids(candidate_rsids: set[str]) -> set[str]:
    present = set()
    with REF_BIM.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and fields[1] in candidate_rsids:
                present.add(fields[1])
    return present


def prepare():
    WORK_DIR.mkdir(exist_ok=True)
    CLUMP_DIR.mkdir(exist_ok=True)
    SCORE_DIR.mkdir(exist_ok=True)

    raw = {}
    candidate_rsids = set()
    summary = []
    for trait in TRAITS:
        rows = list(iter_gwas(trait))
        raw[trait] = rows
        candidate_rsids.update(row[0] for row in rows)
        summary.append([trait, "candidate_p_le_5e-8", len(rows)])

    present = ref_rsids(candidate_rsids)
    assoc_path = WORK_DIR / "assoc_p_le_5e-8.tsv"
    with assoc_path.open("w", newline="") as assoc_handle:
        assoc_writer = csv.writer(assoc_handle, delimiter="\t", lineterminator="\n")
        assoc_writer.writerow(["trait", "SNP", "A1", "A2", "BETA", "P"])
        for trait in TRAITS:
            by_snp = {}
            missing_ref = 0
            for snp, a1, a2, beta, p_value in raw[trait]:
                if snp not in present:
                    missing_ref += 1
                    continue
                if snp not in by_snp or p_value < by_snp[snp][4]:
                    by_snp[snp] = (a1, a2, beta, p_value)

            clump_input = WORK_DIR / f"{trait}.clump_input.tsv"
            with clump_input.open("w", newline="") as clump_handle:
                clump_writer = csv.writer(clump_handle, delimiter="\t", lineterminator="\n")
                clump_writer.writerow(["SNP", "A1", "P"])
                for snp, (a1, a2, beta, p_value) in sorted(by_snp.items(), key=lambda item: item[1][3]):
                    clump_writer.writerow([snp, a1, f"{p_value:.12g}"])
                    assoc_writer.writerow([trait, snp, a1, a2, f"{beta:.16g}", f"{p_value:.12g}"])

            summary.append([trait, "not_in_reference", missing_ref])
            summary.append([trait, "clump_input", len(by_snp)])

    with (WORK_DIR / "prepare_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "metric", "value"])
        writer.writerows(summary)


def run_clump():
    CLUMP_DIR.mkdir(exist_ok=True)
    for trait in TRAITS:
        out_prefix = CLUMP_DIR / trait
        cmd = [
            str(PLINK2),
            "--bfile",
            str(REF_PREFIX),
            "--clump",
            str(WORK_DIR / f"{trait}.clump_input.tsv"),
            "--clump-id-field",
            "SNP",
            "--clump-p-field",
            "P",
            "--clump-a1-field",
            "A1",
            "--clump-p1",
            str(MAX_P),
            "--clump-p2",
            str(MAX_P),
            "--clump-r2",
            "0.1",
            "--clump-kb",
            "250",
            "--out",
            str(out_prefix),
        ]
        subprocess.run(cmd, check=True)


def read_assoc():
    assoc = defaultdict(dict)
    with (WORK_DIR / "assoc_p_le_5e-8.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            assoc[row["trait"]][row["SNP"]] = {
                "a1": row["A1"],
                "a2": row["A2"],
                "beta": float(row["BETA"]),
                "p": float(row["P"]),
            }
    return assoc


def read_clumped(trait: str) -> set[str]:
    path = CLUMP_DIR / f"{trait}.clumps"
    snps = set()
    with path.open() as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if header is None:
                if fields[0] == "#CHROM":
                    fields[0] = "CHROM"
                header = fields
                continue
            snps.add(fields[header.index("ID")])
    return snps


def assemble():
    SCORE_DIR.mkdir(exist_ok=True)
    assoc = read_assoc()
    clumped = {trait: read_clumped(trait) for trait in TRAITS}
    summary = []
    alignment_notes = []

    for label, threshold in THRESHOLDS:
        by_snp = defaultdict(dict)
        for trait in TRAITS:
            for snp in clumped[trait]:
                rec = assoc[trait].get(snp)
                if rec is not None and rec["p"] <= threshold:
                    by_snp[snp][trait] = rec

        output = SCORE_DIR / f"EAS_{label}.score"
        with output.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "A1", *TRAITS])
            for snp in sorted(by_snp):
                records = by_snp[snp]
                _, anchor = min(records.items(), key=lambda item: item[1]["p"])
                score_a1 = anchor["a1"]
                betas = []
                for trait in TRAITS:
                    rec = records.get(trait)
                    if rec is None:
                        betas.append("0.0")
                    elif rec["a1"] == score_a1:
                        betas.append(f"{rec['beta']:.16g}")
                    elif rec["a2"] == score_a1:
                        betas.append(f"{-rec['beta']:.16g}")
                        alignment_notes.append([label, snp, trait, rec["a1"], rec["a2"], score_a1, "beta_flipped"])
                    else:
                        betas.append("0.0")
                        alignment_notes.append([label, snp, trait, rec["a1"], rec["a2"], score_a1, "allele_mismatch_set_zero"])
                writer.writerow([snp, score_a1, *betas])

        counts = {trait: 0 for trait in TRAITS}
        for records in by_snp.values():
            for trait in records:
                counts[trait] += 1
        summary.append([label, "union_snps", len(by_snp)])
        for trait in TRAITS:
            summary.append([label, trait, counts[trait]])

    with (SCORE_DIR / "score_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["threshold", "metric", "value"])
        writer.writerows(summary)

    with (SCORE_DIR / "allele_alignment_notes.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["threshold", "SNP", "trait", "trait_A1", "trait_A2", "score_A1", "action"])
        writer.writerows(alignment_notes)

    for label, _ in THRESHOLDS:
        source = SCORE_DIR / f"EAS_{label}.score"
        target = BASE_DIR / f"EAS_{label}.score"
        target.write_text(source.read_text())


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in {"prepare", "clump", "assemble", "all"}:
        raise SystemExit("Usage: pt_pipeline_eas.py {prepare|clump|assemble|all}")
    if sys.argv[1] in {"prepare", "all"}:
        prepare()
    if sys.argv[1] in {"clump", "all"}:
        run_clump()
    if sys.argv[1] in {"assemble", "all"}:
        assemble()


if __name__ == "__main__":
    main()
