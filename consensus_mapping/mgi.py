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
parser = argparse.ArgumentParser(description='MGI search')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--mgi_data_path", type=str)
parser.add_argument("--gene_annotation_file", type=str)
parser.add_argument("--mpo_dict", type=str)

args = parser.parse_args()

input_path = args.input_path #path of input
output_path = args.output_path #path of output
mgi_data_path = args.mgi_data_path
gene_annotation_file = args.gene_annotation_file #annotation file for genes
mpo_dict = eval(args.mpo_dict) #MGI dictionary
#paramerte--------------------------------------------

#PhenotypicAllele
mgi = pd.read_csv('%s/MGI_PhenotypicAllele.rpt' % (mgi_data_path),sep='\t',header=None,comment='#')
mgi.columns = ['allele_accession_id','allele_symbol','allele_name','allele_type','allele_attribute','pmid',
               'marker_accession_id', 'marker_symbol','marker_refseq','marker_ensembl','top_level_mp_term_id','synonyms',
               'marker_name']
mgi_id_symbol_dict = dict(mgi[['marker_accession_id','marker_symbol']].values)

#GenePheno
mgi_mp = pd.read_csv('%s/MGI_GenePheno.rpt' % (mgi_data_path),sep='\t',header=None,comment='#')
mgi_mp.columns = ['allelic_composition','allele_symbol','allele_accession_id','background','mp_term_id',
                  'pmid','marker_accession_id','genotype_id']
mgi_mp['marker_symbol'] = mgi_mp['marker_accession_id'].apply(lambda x: mgi_id_symbol_dict.get(x))
mgi_mp = mgi_mp.dropna(subset=['marker_symbol'])

#MP ontology
graph = obonet.read_obo('%s/MPheno_OBO.ontology' % (mgi_data_path))
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
graph_with_obs = obonet.read_obo('%s/MPheno_OBO.ontology' % (mgi_data_path), ignore_obsolete=False)
old_to_new = dict()
for node, data in graph_with_obs.nodes(data=True):
    for replaced_by in data.get("replaced_by", []):
        old_to_new[node] = replaced_by
mgi_mp['mp_term_id'] = mgi_mp['mp_term_id'].apply(lambda x: old_to_new.get(x,x))

#gene annotation file
gene = pd.read_csv(gene_annotation_file,sep='\t')


#processing each loci
result = []
loci = pd.read_csv('%s/loci.parameter' % (input_path),sep=' ',header=None)
loci.columns = ['prefix','locus','indsig','chr','s','e']
loci = loci.drop_duplicates(subset=['prefix','locus'])
for i in loci.index:
    prefix,locus,chr_,start,end = loci.loc[i,['prefix','locus','chr','s','e']]
    mp_original = mpo_dict[prefix].split(',')
    #list of eligible MP
    for mp in mp_original:
        mp_trait_lst = [x for x in set(mgi_mp['mp_term_id'].values) if np.any([i in [mp] for i in networkx.descendants(graph, x)])] + [mp]
        #select gene
        gene_temp = gene[(gene['CHR']=='chr'+str(chr_)) & ~((gene['END']<start-1000000) | (gene['START']>end+1000000))]
        gene_temp['Gene_mouse'] = gene_temp['NAME'].apply(lambda x: x[0]+x[1::].lower())
        #merge with MGI GenePheno
        gene_temp = pd.merge(gene_temp[['ENSGID','NAME','Gene_mouse']],mgi_mp,left_on='Gene_mouse',right_on='marker_symbol',
                             how='inner')
        #loop each allele
        for mgi_id in gene_temp['allele_accession_id'].unique():
            mgi_temp = gene_temp[gene_temp['allele_accession_id']==mgi_id]
            ensgid,symbol,mouse_symbol = mgi_temp[['ENSGID','NAME','marker_symbol']].values[0]
            allele = mgi_temp['allele_symbol'].values[0]
            gene_id = mgi_temp['marker_accession_id'].values[0]
            mp_lst = mgi_temp['mp_term_id'].values
            #exposure and outcome MP list
            mp_trait_lst_ = [x for x in mp_lst if x in mp_trait_lst]
            if len(mp_trait_lst_)==0:
                continue
            mp_ = [id_to_name[x] for x in mp_trait_lst_]
            result.append([prefix,locus,ensgid,symbol,mouse_symbol,gene_id,mgi_id,allele,mp,
                           '; '.join(mp_),'; '.join(mp_trait_lst_)])
#save
mgi_df = pd.DataFrame(result,columns=['Prefix','Locus','Gene_id','Gene','Gene_mouse','MGI_gene_id','Allele_id',
                                      'Allele_symbol','MPO_original','MPO','MPO_lst'])
mgi_df.to_csv('%s/MGI.csv' % (output_path),index=False)

































