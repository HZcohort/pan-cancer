# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='Loci verify')

parser.add_argument("--gwas_summary_path", type=str)
parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--prefix", type=str)
parser.add_argument("--mode_pleiotropy", action='store_true')
parser.add_argument("--mode_GWAS", action='store_true')
parser.add_argument("--primo", type=float)
parser.add_argument("--flip_bp", type=int)
args = parser.parse_args()

input_path = args.input_path
gwas_summary_path = args.gwas_summary_path
output_path = args.output_path
prefix = args.prefix #file prefix
mode_pleiotropy = args.mode_pleiotropy
mode_GWAS = args.mode_GWAS
primo = args.primo #file prefix
flip_bp = args.flip_bp
#paramerte--------------------------------------------

#loead merged GWAS
merged = pd.read_csv('%s/%s.gz' % (gwas_summary_path,prefix),sep='\t')
#load information for identified loci
loci = pd.read_csv('%s/%s.bed' % (input_path,prefix),sep='\t')
loci['name'] = loci.apply(lambda row: '%i-%i-%i' % (row[0],row[1],row[2]),axis=1)
# load clumpted results for independent significant SNPs
indsig = pd.read_csv('%s/%s_indsig.clumped' % (input_path,prefix),delim_whitespace=True)

lead_lst = []
lead_failed_lst = []
loci_parameter = []
for i in loci.index:
    chr_,start,end,name = loci.loc[i,['0','1','2','name']]
    #merged GWAS information
    temp = merged[(merged['chr']==chr_) & (merged['bp']>=start) & (merged['bp']<=end)]
    temp = temp.sort_values(by='p',ascending=True)
    temp['Loci'] = name
    #lead SNP
    lead,lead_bp = temp.iloc[0][['snp','bp']]
    temp = temp.iloc[[0]]
    #independent significant SNP information
    indsig_locus = indsig[(indsig['CHR']==chr_) & (indsig['BP']>=start) & (indsig['BP']<=end)]
    # if len(indsig_locus) == 1 and indsig_locus.iloc[0]['SNP'] == lead:
    #     indsig_snp = ''
    #     #temp['indsig_snp'] = ''
    # else:
    #     indsig_snp = indsig_locus['SNP'].to_list()
    #     #temp['indsig_snp'] = ';'.join(indsig_locus['SNP'].to_list())
    
    #select locus with at leat one indsig_PP > cut-off
    if mode_pleiotropy:
        save_cols = ['Loci','snp','indsig_snp_pp','chr','bp','p','PP','ref','alt','N',
                     'A1','A2','beta_hat_1','beta_hat_2','seb1','seb2','p1','p2']
        indsig_snps = indsig_locus['SNP'].to_list()
        indsig_snps_pp = merged[merged['snp'].isin(indsig_snps)]
        temp['indsig_snp_pp'] = [dict(indsig_snps_pp[['snp','PP']].values)]
        if np.any(indsig_snps_pp['PP'] >= primo):
            lead_lst.append(temp[save_cols])
            #save paramters for shell command
            delta = flip_bp
            start -= delta; end += delta
            snp_lst = indsig_snps_pp[indsig_snps_pp['PP']>=primo]['snp'].to_list()
            for snp in snp_lst:
                loci_parameter.append([prefix,name,snp,chr_,start,end])
        else:
            lead_failed_lst.append(temp[save_cols])
    if mode_GWAS:
        save_cols = ['Loci','snp','indsig_snp_p','chr','bp','p','ref','alt','N',
                     'A1','A2','beta','se']
        indsig_snps = indsig_locus['SNP'].to_list()
        indsig_snps_p = merged[merged['snp'].isin(indsig_snps)]
        temp['indsig_snp_p'] = [dict(indsig_snps_p[['snp','p']].values)]
        #save
        lead_lst.append(temp[save_cols])
        #save paramters for shell command
        delta = flip_bp
        start -= delta; end += delta
        for snp in indsig_snps:
            loci_parameter.append([prefix,name,snp,chr_,start,end])

if mode_pleiotropy:
    try:
        lead_df = pd.concat(lead_lst)
    except:#if empty list
        lead_df = pd.DataFrame([],columns=save_cols)
    try:
        lead_failed_df = pd.concat(lead_failed_lst)
    except:#if empty list
        lead_failed_df = pd.DataFrame([],columns=save_cols)
    bed = loci[loci['name'].isin(lead_df['Loci'])]
    lead_df.to_csv('%s/%s_lead.csv' % (output_path,prefix),index=False)
    lead_failed_df.to_csv('%s/%s_failed_lead.csv' % (output_path,prefix),index=False)
    bed[['0','1','2']].sort_values(by='0').to_csv('%s/%s.bed' % (output_path,prefix),
                                                  index=False,sep='\t')
    pd.DataFrame(loci_parameter,columns=None).to_csv('%s/%s_loci.parameter' % (output_path,prefix),sep=' ',
                                                     header=False,index=False)
if mode_GWAS:
    lead_df = pd.concat(lead_lst)
    bed = loci[loci['name'].isin(lead_df['Loci'])]
    lead_df.to_csv('%s/%s_lead.csv' % (output_path,prefix),index=False)
    bed[['0','1','2']].sort_values(by='0').to_csv('%s/%s.bed' % (output_path,prefix),
                                                  index=False,sep='\t')
    pd.DataFrame(loci_parameter,columns=None).to_csv('%s/%s_loci.parameter' % (output_path,prefix),sep=' ',
                                                     header=False,index=False)



































