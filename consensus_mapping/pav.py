# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='PAV generate file')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--bim_file", type=str)
parser.add_argument("--flip_kb", type=int)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
bim_file = args.bim_file #path of common data
flip_kb = args.flip_kb
#paramerte--------------------------------------------

parameter_lst = []
bim_1000g = pd.read_csv(bim_file,sep='\t',header=None)
parameter = pd.read_csv('%s/loci.parameter' % (input_path),sep=' ',header=None)
parameter.columns = ['prefix','locus','indsig','chr','s','e']

for prefix,locus in parameter[['prefix','locus']].value_counts().index:
    temp = parameter[(parameter['prefix']==prefix) & (parameter['locus']==locus)]
    indsig_lst = temp['indsig'].to_list()
    chr_,start,end = temp.iloc[0][['chr','s','e']]
    #independent significant SNP
    indsig_df = pd.DataFrame(indsig_lst,columns=['SNP'])
    indsig_df['P'] = 5e-9
    #bim file
    bim_temp = bim_1000g[(bim_1000g[0]==chr_) & (bim_1000g[3]>=start-(flip_kb*1000)) & (bim_1000g[3]<=end+(flip_kb*1000))]
    bim_temp['SNP'] = bim_temp[1]
    bim_temp['P'] = 5e-9
    #save
    parameter_lst.append([prefix,locus])
    indsig_df.to_csv('%s/%s_%s.indsig' % (output_path,prefix,locus),index=False,sep='\t')
    bim_temp[['SNP','P']].to_csv('%s/%s_%s.sumstats' % (output_path,prefix,locus),index=False,sep='\t')

#save parameter
pd.DataFrame(parameter_lst).to_csv('%s/pav.parameter' % (output_path),index=False,header=False,sep=' ')
























