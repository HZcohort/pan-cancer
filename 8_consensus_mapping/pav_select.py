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


proxy = pd.read_csv('%s/proxy_snp.tsv' % (input_path),sep='\t')
vep = pd.read_csv('%s/VEP_selected.tsv' % (input_path),sep='\t',header=None)
vep.columns = ['proxy_snp','Location','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position',
               'CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','IMPACT','DISTANCE','STRAND',
               'FLAGS','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL','MANE','TSL','APPRIS','ENSP','GENE_PHENO',
               'SIFT','PolyPhen','AF','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS',
               'MOTIF_SCORE_CHANGE','TRANSCRIPTION_FACTORS']
vep['Gene_id'] = vep['Gene']
vep['Gene'] = vep['SYMBOL']
vep_pav = vep[vep['IMPACT'].isin(['MODERATE','HIGH'])]
vep_other = vep[~vep['IMPACT'].isin(['MODERATE','HIGH'])]

#merge and save
save_cols = ['proxy_snp','Allele','Location','Consequence','IMPACT','Gene_id','Gene',
             'BIOTYPE','Codons','CLIN_SIG','SOMATIC','PHENO','PUBMED']
proxy_pav = pd.merge(proxy,vep_pav[save_cols],on='proxy_snp',how='inner')
proxy_other = pd.merge(proxy,vep_other[save_cols],on='proxy_snp',how='inner')
proxy_pav.drop_duplicates().to_csv('%s/PAV_eligible.csv' % (output_path),index=False)
proxy_other.drop_duplicates().to_csv('%s/PAV_other.csv' % (output_path),index=False)





























