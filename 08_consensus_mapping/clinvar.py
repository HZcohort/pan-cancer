# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import networkx
import obonet
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='ClinVar search')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--clinvar_data_path", type=str)
parser.add_argument("--gene_annotation_file", type=str)
parser.add_argument("--hpo_dict", type=str)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
clinvar_data_path = args.clinvar_data_path
gene_annotation_file = args.gene_annotation_file #annotation file for genes
hpo_dict = eval(args.hpo_dict) #HPO dictionary
#paramerte--------------------------------------------


#processing Clinvar DB
clinvar = pd.read_csv('%s/clinvar_GRCh37.vcf.gz' % (clinvar_data_path),sep='\t',comment='#') #read DB
clinvar.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
clinvar['dict'] = clinvar['INFO'].apply(lambda x: {item.split('=')[0]:item.split('=')[1] for item in x.split(';')})
for var in ['ALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG']:
    clinvar[var] = clinvar['dict'].apply(lambda x: x.get(var))
clinvar = clinvar.dropna(subset=['CLNDISDB'])
clinvar = clinvar[['Pathogenic' in x for x in clinvar['CLNSIG']]]
clinvar = clinvar[~clinvar['CLNREVSTAT'].isin(['no_assertion_provided','no_assertion_criteria_provided'])]
#MEDGENE
clinvar['medgen_list'] = clinvar['CLNDISDB'].apply(lambda x: [i.split(':')[1] for s in x.split('|') for i in s.split(',') if 'MedGen' in i])

#processing Network of HPO ontology
graph = obonet.read_obo('%s/hp.obo' % (clinvar_data_path))
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
graph_with_obs = obonet.read_obo('%s/hp.obo' % (clinvar_data_path), ignore_obsolete=False)
old_to_new = dict()
for node, data in graph_with_obs.nodes(data=True):
    for replaced_by in data.get("replaced_by", []):
        old_to_new[node] = replaced_by

#MedGen mapping to HPO
medgen_hpo = pd.read_csv('%s/MedGen_HPO_OMIM_Mapping.txt.gz' % (clinvar_data_path),sep='|')
#medgen_hpo = medgen_hpo[medgen_hpo['#OMIM_CUI'].isin(medgene_list)]
medgen_hpo['HPO_ID'] = medgen_hpo['HPO_ID'].apply(lambda x: old_to_new.get(x,x))
medgen_name_dict = dict(medgen_hpo[['#OMIM_CUI','OMIM_name']].values)
hpo_name_dict = dict(medgen_hpo[['HPO_ID','HPO_name']].values)

#HPO list for each MedGen
medgen_hpo_dict = {}
for medgen in medgen_hpo['#OMIM_CUI'].unique():
    hpo_lst = list(set(medgen_hpo[medgen_hpo['#OMIM_CUI']==medgen]['HPO_ID'].to_list()))
    medgen_hpo_dict[medgen] = hpo_lst

#ALLELEID to gene
Clinvar_gene = pd.read_csv('%s/allele_gene.txt.gz' % (clinvar_data_path),sep='\t')
Clinvar_gene.rename(columns={'#AlleleID':'ALLELEID','Symbol':'Gene'},inplace=True)

#gene annotation file
gene = pd.read_csv(gene_annotation_file,sep='\t')
gene_id_dict = dict(gene[['NAME','ENSGID']].values)

#processing each loci
df_result = []
loci = pd.read_csv('%s/loci.parameter' % (input_path),sep=' ',header=None)
loci.columns = ['prefix','locus','indsig','chr','s','e']
loci = loci.drop_duplicates(subset=['prefix','locus'])
for i in loci.index:
    prefix,locus,chr_,start,end = loci.loc[i,['prefix','locus','chr','s','e']]
    hpo_orginal_lst = hpo_dict[prefix].split(',')
    #dictionary of eligible HPO
    for hpo in hpo_orginal_lst:
        hpo_lst = set([x for x in set(medgen_hpo['HPO_ID'].values) if np.any([i in [hpo] for i in networkx.descendants(graph, x)])] + [hpo])
        #subset of clinvar
        clinvar_df = clinvar[(clinvar['CHROM']==chr_) & (clinvar['POS']>=start-1000000) & (clinvar['POS']<=end+1000000)]
        #get the HPO list for each MedGen
        clinvar_df['HPO_lst'] = clinvar_df['medgen_list'].apply(lambda x: [hpo for medgen in x for hpo in medgen_hpo_dict.get(medgen,[])])
        #HPO list related to the trait
        clinvar_df['HPO_lst'] = clinvar_df['HPO_lst'].apply(lambda x: list(set(x).intersection(hpo_lst)))
        #keep length of HPO_lst list > 0
        clinvar_df = clinvar_df[[len(x)>0 for x in clinvar_df['HPO_lst']]]
        if clinvar_df.empty:
            continue
        #HPO name
        clinvar_df['HPO'] = clinvar_df['HPO_lst'].apply(lambda x: '; '.join([hpo_name_dict[i] for i in x]) if len(x)>0 else np.nan)
        #megen name
        clinvar_df['MedGene_name'] = clinvar_df['medgen_list'].apply(lambda x: [medgen_name_dict.get(i) for i in x])
        clinvar_df['MedGene_name'] = clinvar_df['MedGene_name'].apply(lambda x: '; '.join([i for i in x if i is not None]))
        #mapping to gene
        clinvar_df['ALLELEID'] = clinvar_df['ALLELEID'].astype(int)
        clinvar_df = pd.merge(clinvar_df,Clinvar_gene[['ALLELEID','Gene']],on=['ALLELEID'],how='inner')
        #save
        clinvar_df['Prefix'] = prefix
        clinvar_df['Locus'] = locus
        clinvar_df['Gene_id'] = clinvar_df['Gene'].apply(lambda x: gene_id_dict.get(x))
        clinvar_df['HPO_original'] = hpo #original HPO
        save_cols = ['Prefix', 'Locus','Gene_id','Gene','ALLELEID', 'CHROM','POS','REF','ALT', 'Gene', 'CLNSIG', 'CLNREVSTAT', 
                     'HPO_original','MedGene_name', 'HPO', 'medgen_list', 'HPO_lst']
        df_result.append(clinvar_df[save_cols])

pd.concat(df_result).to_csv('%s/ClinVar.csv' % (output_path),index=False)





































