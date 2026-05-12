# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

#paramerte----------------------
parser = argparse.ArgumentParser(description='Identiify the loci')

parser.add_argument("--input_path", type=str)
parser.add_argument("--g1000_ref", type=str)
parser.add_argument("--regions_merge_bp", type=int)
args = parser.parse_args()

input_path = args.input_path #path of project
g1000_ref = args.g1000_ref
regions_merge_bp = args.regions_merge_bp
#paramerte----------------------

bim_1000g = pd.read_csv(g1000_ref,sep='\t',header=None)
bp_dict = dict(bim_1000g[[1,3]].values)

#generate bed files
for dir_,_,files in os.walk(input_path):
    for file in files:
        if 'locus.clumped' in file:
            path = os.path.join(dir_,file)
            prefix = file.split('.')[0].replace('_locus','')
            print(prefix)
            clump = pd.read_csv(path,delim_whitespace=True)
            clump['rs'] = clump['SP2'].apply(lambda x: [s.split('(')[0] for s in x.split(',')])
            clump['min_bp'] = clump['rs'].apply(lambda x: min([bp_dict[s] for s in x]))
            clump['max_bp'] = clump['rs'].apply(lambda x: max([bp_dict[s] for s in x]))
            bed = []
            for chr_ in clump['CHR'].unique():
                temp = clump[clump['CHR']==chr_]
                if len(temp) == 1:
                    bed.append([chr_,temp['min_bp'].values[0],temp['max_bp'].values[0]])
                else:
                    temp = temp.sort_values(by=['min_bp','max_bp'])
                    c_min,c_max = 0,0
                    for i in temp.index:
                        min_ = temp.loc[i,'min_bp']
                        max_ = temp.loc[i,'max_bp']
                        if c_min+c_max == 0:
                            c_min, c_max = min_, max_
                            continue
                        elif min_ <= c_max + regions_merge_bp:
                            c_max =  max(max_,c_max)
                            continue
                        elif min_ > c_max + regions_merge_bp:
                            bed.append([chr_,c_min,c_max])
                            c_min,c_max = min_,max_
                            continue
                    if c_min+c_max > 0:
                        bed.append([chr_,c_min,c_max])
            pd.DataFrame(bed).sort_values(by=0).to_csv('%s/%s.bed' % (input_path,prefix),index=False,sep='\t')

















































