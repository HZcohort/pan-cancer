# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='PoPS merge')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--gene_annot_file", type=str)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
gene_annot_file = args.gene_annot_file
#paramerte--------------------------------------------

parameter_file = pd.read_csv('%s/MAGMA/PoPS.parameter' % (input_path),sep=' ',header=None)
parameter_file.columns = ['prefix','loci','n']
gene_file = pd.read_csv(gene_annot_file,sep='\t')
id_symbol_dict = dict(gene_file.dropna()[['ENSGID','NAME']].values)

df_result = []
save_cols = ['Prefix','Locus','Gene_id','Gene','PoPS_Score','Rank','Feature1','Feature2','Feature3','Feature4','Feature5',
             'Feature6','Feature7','Feature8','Feature9','Feature10']
for i in parameter_file.index:
    prefix,loci = parameter_file.loc[i,['prefix','loci']]
    pops = pd.read_csv('%s/data/%s.preds' % (input_path,prefix),sep='\t').sort_values(by=['PoPS_Score'],ascending=False)
    pops = pd.merge(pops,gene_file[['ENSGID','NAME','CHR','START','END']],on='ENSGID',how='left')
    feature = np.load('%s/data/%s_feature.npy' % (input_path,prefix))
    for locus in loci.split(','):
        chr_,s,e = [int(x) for x in locus.split('-')]
        pops_temp = pops[(pops['CHR']==chr_) & (pops['START']<=e+1000000) & (pops['END']>=s-1000000)]
        try:
            pops_temp = pops_temp[['ENSGID','NAME','PoPS_Score']]
        except:
            pops_temp = pops_temp[['ENSGID','NAME','PoPS_Score']]
        feature_index = np.where([x in pops_temp.ENSGID.values for x in feature[:,0]])[0]
        feature_ = pd.DataFrame(feature[feature_index,0:11],columns=['ENSGID']+['Feature%i' % i for i in range(1,11)])
        pops_temp = pd.merge(pops_temp,feature_,on='ENSGID',how='left')
        pops_temp['Prefix'] = prefix
        pops_temp['Locus'] = locus
        pops_temp['Rank'] = np.arange(1,len(pops_temp)+1)
        pops_temp.rename(columns={'ENSGID':'Gene_id','NAME':'Gene'},inplace=True)
        df_result.append(pops_temp[save_cols])

#save
pd.concat(df_result).to_csv('%s/PoPS.csv' % (output_path),index=False)
























