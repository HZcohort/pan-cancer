# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='PAV clump')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
#paramerte--------------------------------------------

result = []
parameter = pd.read_csv('%s/pav.parameter' % (input_path),header=None,sep=' ')
for i in parameter.index:
    prefix,locus = parameter.loc[i,[0,1]]
    clump = pd.read_csv('%s/%s_%s.clumped' % (input_path,prefix,locus),delim_whitespace=True)
    clump['rs'] = clump['SP2'].apply(lambda x: [s.split('(')[0] for s in x.split(',')])
    for j in clump.index:
        indsig = clump.loc[j,'SNP']
        rs_lst = clump.loc[j,'rs']
        for rs in rs_lst:
            result.append([prefix,locus,indsig,rs])
#file for proxy snp
df_result = pd.DataFrame(result,columns=['prefix','loci','indsig_snp','proxy_snp'])
#save
df_result[['proxy_snp']].to_csv('%s/proxy_snp.list' % (output_path),index=False,header=False,sep='\t')
df_result.to_csv('%s/proxy_snp.tsv' % (output_path),index=False,sep='\t')












