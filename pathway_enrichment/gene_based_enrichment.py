# -*- coding: utf-8 -*-

import argparse
import json
from pathlib import Path

import gseapy as gp
import pandas as pd


def read_gene_set(path):
    with open(path, encoding="utf-8") as f:
        gene_set = json.load(f)
    return {pathway: [gene for gene in gene_set[pathway]["geneSymbols"]] for pathway in gene_set.keys()}


parser = argparse.ArgumentParser(description="Run gene-based enrichment for mapped pan-cancer risk genes")
parser.add_argument("--gene_id_map", required=True, help="TSV with a NAME column used as the background gene universe")
parser.add_argument("--msigdb_json_path", required=True, help="Directory containing MSigDB JSON files")
parser.add_argument("--magma_gene_set", required=True, help="Merged MAGMA gene-set output from MAGMA_geneset_merge.py")
parser.add_argument("--consensus_genes", required=True, help="Consensus mapped-gene table with a Gene column")
parser.add_argument("--output", required=True, help="Output CSV path")
parser.add_argument("--q_threshold", type=float, default=0.05, help="MAGMA gene-set q-value threshold")
args = parser.parse_args()

msigdb_json_path = Path(args.msigdb_json_path)

gene_lst = pd.read_csv(args.gene_id_map, sep="\t")

go = read_gene_set(msigdb_json_path / "c5.go.v7.5.1.Hs.json")
cp = read_gene_set(msigdb_json_path / "c2.cp.v7.5.1.Hs.json")
hallmark = read_gene_set(msigdb_json_path / "h.all.v7.5.1.Hs.json")

sig_path = pd.read_csv(args.magma_gene_set)
sig_path = sig_path[sig_path["Q"] <= args.q_threshold]
sig_cp = sig_path[sig_path["source"].isin(["C2:CP:BIOCARTA", "C2:CP", "C2:CP:REACTOME", "C2:CP:WIKIPATHWAYS", "C2:CP:PID"])]
sig_gobp = sig_path[sig_path["source"].isin(["C5:GO:BP"])]
sig_h = sig_path[sig_path["source"].isin(["H"])]

gobp_dict = {i: j for i, j in go.items() if i in sig_gobp["FULL_NAME"].values}
BIOCARTA_dict = {i: j for i, j in cp.items() if i in sig_cp[sig_cp["source"].isin(["C2:CP:BIOCARTA"])]["FULL_NAME"].values}
CP_dict = {i: j for i, j in cp.items() if i in sig_cp[sig_cp["source"].isin(["C2:CP"])]["FULL_NAME"].values}
REACTOME_dict = {i: j for i, j in cp.items() if i in sig_cp[sig_cp["source"].isin(["C2:CP:REACTOME"])]["FULL_NAME"].values}
WIKIPATHWAYS_dict = {i: j for i, j in cp.items() if i in sig_cp[sig_cp["source"].isin(["C2:CP:WIKIPATHWAYS"])]["FULL_NAME"].values}
PID_dict = {i: j for i, j in cp.items() if i in sig_cp[sig_cp["source"].isin(["C2:CP:PID"])]["FULL_NAME"].values}
H_dict = {i: j for i, j in hallmark.items() if i in sig_h[sig_h["source"].isin(["H"])]["FULL_NAME"].values}

consensus = pd.read_csv(args.consensus_genes)
gene_mapped_lst = consensus["Gene"].unique().tolist()

enr2 = gp.enrich(
    gene_list=gene_mapped_lst,
    gene_sets=[gobp_dict, BIOCARTA_dict, CP_dict, REACTOME_dict, WIKIPATHWAYS_dict, PID_dict, H_dict],
    background=gene_lst["NAME"].values,
    outdir=None,
    verbose=True,
)

df = enr2.results
gene_set_name_dict = {
    "gs_ind_0": "GO BP",
    "gs_ind_1": "BIOCARTA",
    "gs_ind_2": "Canonical pathways",
    "gs_ind_3": "REACTOME",
    "gs_ind_4": "WIKIPATHWAYS",
    "gs_ind_5": "PID",
    "gs_ind_6": "HALLMARK",
}
df["Gene_set"] = df["Gene_set"].apply(lambda x: gene_set_name_dict.get(x))
total_n = len(df)
df["Adjusted P-value"] = df["P-value"].apply(lambda x: x * total_n if x * total_n <= 1 else 1)
df["n_gene"] = df["Overlap"].apply(lambda x: int(x.split("/")[0]))
df = df.sort_values(by="Adjusted P-value")

df.to_csv(args.output, index=False)
