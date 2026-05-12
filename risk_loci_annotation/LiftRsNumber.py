# -*- coding: utf-8 -*-
#lift the rs in GWAS summary
import os
import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse
import warnings
warnings.filterwarnings("ignore")


#paramerte----------------------
parser = argparse.ArgumentParser(description='Liftover')

parser.add_argument("--input_path", type=str)
parser.add_argument("--trait", type=str)
parser.add_argument("--liftover_file", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--bim_file", type=str)
parser.add_argument("--columns", type=str)
args = parser.parse_args()

input_path = args.input_path 
liftover_file = args.liftover_file 
trait = args.trait
bim_file = args.bim_file
columns = eval(args.columns)
output_path = args.output_path
#paramerte----------------------

"""
def myopen(fn):
    import gzip
    try:
        h = gzip.open(fn)
        ln = h.read(2) # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn)
    h.close()
    return gzip.open(fn)

RS_HISTORY = set() # store rs
RS_MERGE = dict() # high_rs -> (lower_rs, current_rs)
	
for ln in myopen(os.path.join(liftover_file,'SNPHistory.bcp.gz')):
    ln = ln.decode()
    fd = ln.split('\t')
    if ln.find('re-activ') < 0:
        RS_HISTORY.add(fd[0])

# record rs number merge history
for ln in myopen(os.path.join(liftover_file,'RsMergeArch.bcp.gz')):
	fd = ln.decode().split('\t')
	h, l = fd[0], fd[1]
	c = fd[6]
	# print 'c=', c
	RS_MERGE[h] = (l, c)

np.save(os.path.join(liftover_file,'rs_history.npy'),RS_HISTORY)
np.save(os.path.join(liftover_file,'rs_merge.npy'),RS_MERGE)
"""

"""
#select specific trait
file_excel = pd.read_excel(os.path.join(project_path,'result/parameter_LDSC.xlsx'))
file_excel = file_excel[file_excel['Traits']==trait]


#loop over all the traits to check information first
for i in file_excel.index:
    trait,snp,a1,a2,p,beta,OR,se,n = file_excel.loc[i,['Traits','snp', 'a1', 'a2', 'p','beta', 'or', 'se', 'N']]
    if np.any(pd.isna([trait,snp,a1,a2,n])):
        raise Exception('Missing information for trait line %i' % (i+1))
    elif pd.isna(beta) and pd.isna(OR):
        raise Exception('Missing information for trait %s: effect size' % (trait))
    elif pd.isna(p) and pd.isna(se):
        raise Exception('Missing information for trait %s: se and p-values' % (trait))
    else:
        continue
print('Check complete for trait %s' % (trait))
"""

#read GWAS

snp,a1,a2,p,beta,OR,se,n = columns #column name
df = pd.read_csv(os.path.join(input_path,'%s.gz' % (trait)),sep='\t')
#first deal with OR
if beta=='':
    df[OR] = df[OR].apply(pd.to_numeric, errors='coerce')
    print('%s: Exclude invalid OR: %i' % (trait,len(df[(df[OR]<0) | (df[OR].isna()) | (df[OR]==1)])))
    df = df[(df[OR]>0) & (~df[OR].isna()) & (df[OR]!=1)]
    df['beta'] = np.log(df[OR])
    beta = 'beta'
else:
    df[beta] = df[beta].apply(pd.to_numeric, errors='coerce')
    print('%s: Exclude invalid beta: %i' % (trait,len(df[(df[beta]==0) | (df[beta].isna())])))
    df = df[(df[beta]!=0) & (~df[beta].isna())]
# se or p-value is na
if se=='':
    df[p] = df[p].apply(pd.to_numeric, errors='coerce')
    print('%s: Exclude invalid p-value: %i' % (trait,len(df[(df[p]<0) | (df[p]>=1) | (df[p].isna())])))
    df = df[(df[p]>0) & (df[p]<1) & (~df[p].isna())]
    df['Z'] = norm.ppf(df[p]/2)
    df['se'] = np.abs(df[beta]/df['Z'])
    se = 'se'
elif p=='':
    df[se] = df[se].apply(pd.to_numeric, errors='coerce')
    print('%s: Exclude invalid se: %i' % (trait,len(df[(df[se].isna()) | (df[se]<=0)])))
    df = df[(~df[se].isna()) & (df[se]>0)]
    df['Z'] = df[beta]/df[se]
    df['p'] = norm.sf(np.abs(df['Z']))*2
    p = 'p'
else:
    df[p] = df[p].apply(pd.to_numeric, errors='coerce')
    df[se] = df[se].apply(pd.to_numeric, errors='coerce')
    print('%s: Exclude invalid se or p: %i' % (trait,len(df[(df[se].isna()) | (df[se]<=0) | (df[p]<0) | (df[p]>1) | (df[p].isna())])))
    df = df[(~df[se].isna()) & (df[se]>0)]
    df = df[(df[p]>0) & (df[p]<=1) & (~df[p].isna())]
#update the dictionary
col_dict = {snp:'snp',a1:'a1',a2:'a2',p:'p',beta:'beta',se:'se'}
df.rename(columns = col_dict,inplace=True)
df_standard = df[['snp','a1','a2','p','beta','se']]
del df

#liftover
RS_HISTORY = np.load(os.path.join(liftover_file,'rs_history.npy'),allow_pickle=True).item() # store rs
RS_MERGE = np.load(os.path.join(liftover_file,'rs_merge.npy'),allow_pickle=True).item() # high_rs -> (lower_rs, current_rs)
map_dict = {}
rs_lst = df_standard['snp'].astype('str').apply(lambda x: x.replace('rs',''))
for rs in rs_lst:
    try:
        rs = rs.replace('rs','')
    except:
        continue
	#if not rsPattern.match(rs):
	    #print ('ERROR: rs number should be like "1000"')
	    #sys.exit(2)
	# rs number not appear in RS_MERGE -> there is no merge on this rs
    if rs not in RS_MERGE:
	    continue
	# lift rs number
    while True:
        if rs in RS_MERGE:
            rsLow, rsCurrent = RS_MERGE[rs]
            if rsCurrent not in RS_HISTORY:
                map_dict['rs'+rs] = 'rs'+rsCurrent
                break
            else:
                rs = rsLow
        else:
            break
print('Number of variants lifted for trait %s: %i' % (trait,len(map_dict)))
df_standard['snp'] = df_standard['snp'].apply(lambda x: map_dict.get(x,x))
df_standard.sort_values(by='p',inplace=True)
df_standard = df_standard.drop_duplicates(subset=['snp'],keep='first')
#delete to save memory
del RS_HISTORY
del RS_MERGE
del rs_lst

#filter for 1kg snps
bim_kg = pd.read_csv(bim_file,sep='\t',header=None,usecols=[0,1,3,4,5])
bim_kg = bim_kg[~((bim_kg[0]==6) & (bim_kg[3]>=26000000) & (bim_kg[3]<=34000000))]
df_merged = pd.merge(df_standard,bim_kg,left_on='snp',right_on=1,how='inner')
#delete to save memory
del bim_kg
del df_standard
#filter snps consistent in both files
df_merged['a1'] = df_merged['a1'].apply(lambda x: x.upper() if isinstance(x, str) else np.nan)
df_merged['a2'] = df_merged['a2'].apply(lambda x: x.upper() if isinstance(x, str) else np.nan)

def flip_vectorized(series):
    return series.replace({'T': 'A', 'A': 'T', 'C': 'G', 'G': 'C'})

cond1 = (df_merged['a1'] == df_merged[4]) & (df_merged['a2'] == df_merged[5])
cond2 = (df_merged['a1'] == df_merged[5]) & (df_merged['a2'] == df_merged[4])
cond3 = (df_merged['a1'] == flip_vectorized(df_merged[4])) & (df_merged['a2'] == flip_vectorized(df_merged[5]))
cond4 = (df_merged['a1'] == flip_vectorized(df_merged[5])) & (df_merged['a2'] == flip_vectorized(df_merged[4]))

df_merged = df_merged.loc[(cond1 | cond2 | cond3 | cond4)]
del cond1
del cond2
del cond3
del cond4

#save
df_merged.rename(columns = {0:'chr',3:'bp',4:'alt',5:'ref','a1':'A1','a2':'A2'},inplace=True)
df_merged['N'] = int(n)
save_col = ["snp","chr","bp","ref","alt","N","beta","se","A1","A2","p"]
df_merged[save_col].dropna().to_csv(os.path.join(output_path,'%s.gz' % (trait)),sep='\t',index=False)




