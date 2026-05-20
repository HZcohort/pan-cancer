# -*- coding: utf-8 -*-
"""
Merge LDSC genetic-correlation results from phenome-wide analyses.
"""

import os
import pandas as pd
import numpy as np
import argparse
import re
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description="Merge LDSC genetic-correlation results")
parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--prefix", type=str, default=None)
args = parser.parse_args()

input_path = args.input_path
output_path = args.output_path
prefix = args.prefix


def gc_from_ldsc_text(ldsc_file):
    ldsc_result_text = ""
    save_flag = False
    with open(ldsc_file) as f:
        for line in f.readlines():
            if line == "Summary of Genetic Correlation Results\n":
                save_flag = True
                continue
            elif "Analysis finished at" in line:
                save_flag = False
                continue
            if save_flag == True:
                ldsc_result_text += line
            else:
                continue
    ldsc_result_df = []
    ldsc = ldsc_result_text.split("\n")
    for line in ldsc:
        if line == "":
            continue
        else:
            ldsc_result_df.append(line.split())
    ldsc_df = pd.DataFrame(ldsc_result_df[1::], columns=ldsc_result_df[0])
    ldsc_df["p1"] = ldsc_df["p1"].apply(lambda x: x.split("/")[-1].split(".")[0])
    ldsc_df["p2"] = ldsc_df["p2"].apply(lambda x: x.split("/")[-1].split(".")[0])
    return ldsc_df


def h2_lambda_from_ldsc_text(ldsc_file):
    pattern_intercept = r"Intercept: (\d+\.\d*) \((\d+\.\d+)\)"
    pattern_h2_gc = r"Total Observed scale h2: (-?\d*\.?\d+(?:e[+-]?\d+)?) \((\d*\.?\d+(?:e[+-]?\d+)?)\)\nLambda GC: (\d*\.?\d+(?:e[+-]?\d+)?)"
    pattern_heri = r"Heritability of phenotype (\d+)"
    with open(ldsc_file) as f:
        ldsc_text = f.readlines()
    text_dict = {}
    number = None
    for line in ldsc_text:
        if len(re.findall(pattern_heri, line)) > 0:
            number = int(re.findall(pattern_heri, line)[0]) - 1
        if number is not None:
            try:
                text_dict[number] += line
            except:
                text_dict[number] = line
    parameter_dict = {}
    for key in text_dict.keys():
        text = text_dict[key]
        matches_intercept = re.findall(pattern_intercept, text)[0]
        matches_h2_gc = re.findall(pattern_h2_gc, text)[0]
        parameter_dict[key] = [float(x) for x in matches_h2_gc + matches_intercept]
    h2_lambda_intercept_df = pd.DataFrame([], columns=["h2", "h2_se", "lambda", "intercept", "intercept_se"])
    for key in parameter_dict.keys():
        h2_lambda_intercept_df.loc[key] = parameter_dict[key]
    return h2_lambda_intercept_df


def merge_ldsc_h2_lambda_df(ldsc_df, h2_lambda_intercept_df):
    for col in ["h2", "h2_se", "lambda", "intercept", "intercept_se"]:
        try:
            ldsc_df["p1_" + col] = h2_lambda_intercept_df.loc[0, col]
        except:
            ldsc_df["p1_" + col] = np.NaN
    for i in ldsc_df.index:
        index = i + 1
        for col in ["h2", "h2_se", "lambda", "intercept", "intercept_se"]:
            try:
                ldsc_df.loc[i, "p2_" + col] = h2_lambda_intercept_df.loc[index, col]
            except:
                ldsc_df.loc[i, "p2_" + col] = np.NaN
    ldsc_df = ldsc_df.rename(columns={"rg": "rg_ldsc", "se": "rg_se_ldsc", "p": "rg_p_ldsc",
                                      "p2_h2": "p2_h2_ldsc", "p2_h2_se": "p2_h2_se_ldsc",
                                      "p1_h2": "p1_h2_ldsc", "p1_h2_se": "p1_h2_se_ldsc"})
    cols_save = ["p1", "p2", "rg_ldsc", "rg_se_ldsc", "rg_p_ldsc", "gcov_int", "gcov_int_se",
                 "p1_h2_ldsc", "p1_h2_se_ldsc", "p1_lambda", "p1_intercept", "p1_intercept_se",
                 "p2_h2_ldsc", "p2_h2_se_ldsc", "p2_lambda", "p2_intercept", "p2_intercept_se"]
    return ldsc_df[cols_save]


ldsc_lst = []
for dir_, _, files in os.walk(input_path):
    for file in files:
        if "ldsc" in file:
            ldsc_path = os.path.join(dir_, file)
            ldsc_df = gc_from_ldsc_text(ldsc_path)
            h2_lambda_intercept_df = h2_lambda_from_ldsc_text(ldsc_path)
            ldsc_merged = merge_ldsc_h2_lambda_df(ldsc_df, h2_lambda_intercept_df)
            ldsc_lst.append(ldsc_merged)

ldsc_result_df = pd.concat(ldsc_lst)
if prefix:
    ldsc_result_df.to_csv("%s/%s_gc.csv" % (output_path, prefix))
else:
    ldsc_result_df.to_csv("%s/gc.csv" % output_path)
