# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")


#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='H-MAGMA merge')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--parameter_file", type=str)
parser.add_argument("--gene_annot_file", type=str)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
parameter_file = args.parameter_file #H-MAGMA parameters
gene_annot_file = args.gene_annot_file
#paramerte--------------------------------------------

parameter_file = pd.read_csv(parameter_file,sep=' ',header=None)
parameter_file.columns = ['prefix','tissue','loci_lst','n']
gene_file = pd.read_csv(gene_annot_file,sep='\t')
id_symbol_dict = dict(gene_file.dropna()[['ID','Gene']].values)

df_result = []
df_result_sig = []
save_cols = ['Prefix','Locus','Tissue','Gene','Gene_id','START','STOP','NSNPS','ZSTAT','P']
for i in parameter_file.index:
    prefix,tissue,loci_lst = parameter_file.loc[i,['prefix','tissue','loci_lst']]
    out_file = pd.read_csv('%s/%s_%s.genes.out' % (input_path,prefix,tissue),delim_whitespace=True)
    out_file['Gene_id'] = out_file['GENE']
    out_file['Gene'] = out_file['Gene_id'].apply(lambda x: id_symbol_dict.get(x,x))
    #calculate Q-value
    #out_file_sig = out_file[out_file['P']<0.05/len(out_file)]
    #print(i,len(out_file))
    #for each locus
    for locus in loci_lst.split(','):
        chr_,s,e = [int(x) for x in locus.split('-')]
        out_file_temp = out_file[(out_file['CHR']==chr_) & (out_file['START']<=e+1000000) & (out_file['STOP']>=s-1000000)]
        if len(out_file_temp)==0:
            continue
        else:
            out_file_temp['Prefix'] = prefix
            out_file_temp['Locus'] = locus
            out_file_temp['Tissue'] = tissue
            df_result.append(out_file_temp[save_cols])
            #significant result
            out_file_temp_sig = out_file_temp[out_file_temp['P']<0.05/len(out_file)]
            if len(out_file_temp_sig)==0:
                continue
            else:
                df_result_sig.append(out_file_temp_sig[save_cols])

#save all result
try:
    pd.concat(df_result).to_csv('%s/H-MAGMA.csv' % (output_path),index=False)
except:
    pd.DataFrame([],columns=save_cols).to_csv('%s/H-MAGMA.csv' % (output_path),index=False)

#save significant result
try:
    pd.concat(df_result_sig).to_csv('%s/H-MAGMA_significant.csv' % (output_path),index=False)
except:
    pd.DataFrame([],columns=save_cols).to_csv('%s/H-MAGMA_significant.csv' % (output_path),index=False)
















