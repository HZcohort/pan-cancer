# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import argparse
import json
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description="Merge MAGMA gene-set analysis")

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--trait_list", type=str)
parser.add_argument("--geneset_path", type=str)
args = parser.parse_args()

input_path = args.input_path
output_path = args.output_path
trait_lst = args.trait_list
geneset_path = args.geneset_path


def _ecdf(x):
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)


def fdrcorrection(pvals, alpha=0.05):
    """Benjamini-Hochberg FDR correction, adapted from GOATOOLS."""
    pvals = np.asarray(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)

    ecdffactor = _ecdf(pvals_sorted)
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    del pvals_corrected
    return pvals_corrected_


def multiple_testing_correction(ps, alpha=0.05):
    """Correct p-values and return Benjamini-Hochberg q-values."""
    _p = np.array(ps)
    q = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    _q = fdrcorrection(p, alpha)
    q[mask] = _q
    return q


def q_cal(df, prefix_lst):
    q_dict = {}
    for prefix in prefix_lst:
        subset = df[[x.startswith(prefix) for x in df["FULL_NAME"]]]
        subset["Q value"] = multiple_testing_correction(subset["P"].values)
        q_dict.update(dict(subset[["FULL_NAME", "Q value"]].values))
    return q_dict


with open("%s/c2.cp.v7.5.1.Hs.json" % (geneset_path), encoding="utf-8") as f:
    gsea_cp = json.load(f)
with open("%s/c2.cgp.v7.5.1.Hs.json" % (geneset_path), encoding="utf-8") as f:
    gsea_cgp = json.load(f)
with open("%s/c5.go.v7.5.1.Hs.json" % (geneset_path), encoding="utf-8") as f:
    gsea_go = json.load(f)
with open("%s/h.all.v7.5.1.Hs.json" % (geneset_path), encoding="utf-8") as f:
    gsea_hallmark = json.load(f)

info = {}
for dict_ in [gsea_cp, gsea_cgp, gsea_go, gsea_hallmark]:
    info.update({pathway: dict_[pathway] for pathway in dict_.keys()})


sig_lst = []
for trait in trait_lst.split(","):
    sig_lst_trait = []

    cp = pd.read_csv("%s/%s_cp.gsa.out" % (input_path, trait), delim_whitespace=True, comment="#")
    dict_ = q_cal(cp, ["REACTOME", "OTHER", "WP", "KEGG", "BIOCARTA", "PID"])
    cp["Q"] = cp["FULL_NAME"].apply(lambda x: dict_[x])
    cp = cp[["FULL_NAME", "BETA", "BETA_STD", "P", "Q"]]
    sig_lst_trait.append(cp)

    cgp = pd.read_csv("%s//%s_cgp.gsa.out" % (input_path, trait), delim_whitespace=True, comment="#")
    cgp["Q"] = multiple_testing_correction(cgp["P"].values)
    cgp = cgp[["FULL_NAME", "BETA", "BETA_STD", "P", "Q"]]
    sig_lst_trait.append(cgp)

    hallmark = pd.read_csv("%s/%s_hallmark.gsa.out" % (input_path, trait), delim_whitespace=True, comment="#")
    hallmark["Q"] = multiple_testing_correction(hallmark["P"].values)
    hallmark = hallmark[["FULL_NAME", "BETA", "BETA_STD", "P", "Q"]]
    sig_lst_trait.append(hallmark)

    go = pd.read_csv("%s/%s_GO.gsa.out" % (input_path, trait), delim_whitespace=True, comment="#")
    dict_ = q_cal(go, ["GOBP", "GOCC", "GOMF"])
    go["Q"] = go["FULL_NAME"].apply(lambda x: dict_[x])
    go = go[["FULL_NAME", "BETA", "BETA_STD", "P", "Q"]]
    sig_lst_trait.append(go)

    sig = pd.concat(sig_lst_trait)
    sig["msigdb_id"] = sig["FULL_NAME"].apply(lambda x: info[x]["systematicName"])
    sig["source_id"] = sig["FULL_NAME"].apply(lambda x: info[x]["exactSource"])
    sig["url"] = sig["FULL_NAME"].apply(lambda x: info[x]["msigdbURL"])
    sig["source"] = sig["FULL_NAME"].apply(lambda x: info[x]["collection"])
    sig["Prefix"] = trait
    sig_lst.append(sig)

sig_df = pd.concat(sig_lst)[["Prefix", "FULL_NAME", "source", "BETA", "BETA_STD", "P", "Q", "msigdb_id", "source_id", "url"]]
sig_df.to_csv("%s/gene_set_MAGMA.csv" % (output_path))
