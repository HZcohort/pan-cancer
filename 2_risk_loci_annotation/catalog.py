# -*- coding: utf-8 -*-
import pandas as pd
import gzip
import numpy as np
import argparse
from pyliftover import LiftOver
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='Catalog search')

parser.add_argument("--input_path", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--catalog_file", type=str)
parser.add_argument("--liftover_path", type=str)
parser.add_argument("--search_window_kb", type=int, default=100)

args = parser.parse_args()

input_path = args.input_path #path of GWAS file
output_path = args.output_path #path of GWAS file
catalog_file = args.catalog_file #catalog results file
liftover_path = args.liftover_path #path for liftover files
search_window = args.search_window_kb * 1000 #search window (bp) for GWAS catalog database
#paramerte--------------------------------------------

def flip(x):
    if x=='T':
        return 'A'
    elif x=='A':
        return 'T'
    elif x=='C':
        return 'G'
    elif x=='G':
        return 'C'

def bp_pos(x):
    if pd.isna(x):
        return np.NaN
    else:
        try:
            return int(x)
        except:
            try:
                return [int(x) for x in x.split(';')]
            except:
                np.NaN
                
lo = LiftOver('hg19', 'hg38', search_dir=liftover_path) #for liftover

#dealring with catalog file
catalog = pd.read_csv(catalog_file,sep='\t')
catalog = catalog[catalog['P-VALUE']<5e-8]
catalog['CHR_POS'] = catalog['CHR_POS'].apply(lambda x: bp_pos(x))
catalog['CHR_ID'] = catalog['CHR_ID'].apply(lambda x: bp_pos(x))
catalog['index'] = catalog.index
#single snp file
catalog_single = catalog[[type(x) is int for x in catalog['CHR_ID'].values]]
#multi snp
catalog_mult = catalog[[type(x) is list for x in catalog['CHR_ID'].values]]
key_column_lst = [[row[0],row[1][i],row[2][i]] for row in catalog_mult[['index','CHR_ID','CHR_POS']].values 
                  for i in range(len(row[1]))]
key_column_lst = pd.DataFrame(key_column_lst,columns=['index','CHR_ID','CHR_POS'])
catalog_mult_ = pd.merge(catalog_mult[[x for x in catalog_mult.columns if x not in ['CHR_ID','CHR_POS']]],
                         key_column_lst,on='index',how='outer')
catalog_split = pd.concat([catalog_single,catalog_mult_])

#read locus file and search
df_lst = []
loci = pd.read_csv('%s/loci.parameter' % (input_path),sep=' ',header=None)
loci.columns = ['prefix','locus','indsig','chr','s','e']
loci = loci.drop_duplicates(subset=['prefix','locus'])
loci_df_lst = []
for i in loci.index:
    prefix,locus,chr_,start,end = loci.loc[i,['prefix','locus','chr','s','e']]
    #merge loci files
    loci_df = pd.read_csv('%s/%s_lead.csv' % (input_path,prefix))
    cols = ['Prefix'] + list(loci_df.columns)
    loci_df['Prefix'] = prefix
    loci_df_lst.append(loci_df[cols])
    #merge loci files
    try:
        hg38_s = lo.convert_coordinate('chr%i' % chr_,start)[0][1] - search_window
    except: #if cannnot convert, try next cordinates within range
        for bp in range(start,end):
            try:
                hg38_s = lo.convert_coordinate('chr%i' % chr_,bp)[0][1] - search_window
                break
            except:
                continue
    try:
        hg38_e = lo.convert_coordinate('chr%i' % chr_,end)[0][1] + search_window
    except: #if cannnot convert, try next cordinates within range
        for bp in range(end,start,-1):
            try:
                hg38_e = lo.convert_coordinate('chr%i' % chr_,bp)[0][1] + search_window
                break
            except:
                continue
    
    catalog_df = catalog_split[((catalog_split['CHR_ID']==chr_) & (catalog_split['CHR_POS']>=hg38_s) & 
                                (catalog_split['CHR_POS']<=hg38_e))]
    if len(catalog_df) > 0:
        catalog_df['prefix'] = prefix;catalog_df['locus'] = locus
        df_lst.append(catalog_df)
        
catalog_df = pd.concat(df_lst).drop_duplicates(subset=['prefix','locus','index'])
catalog_df = pd.merge(catalog_df[['prefix','locus','index']],catalog,on='index',how='left')
catalog_df.to_csv('%s/catalog.csv' % (output_path),index=False)

loci_df_merged = pd.concat(loci_df_lst).drop_duplicates()
loci_df_merged.to_csv('%s/loci_verified.csv' % (output_path),index=False)




























