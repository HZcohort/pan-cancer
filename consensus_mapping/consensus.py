# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import re
import warnings
warnings.filterwarnings("ignore")

#paramerte--------------------------------------------
parser = argparse.ArgumentParser(description='Consensus based gene mapping')

parser.add_argument("--mode_pleiotropy", action='store_true')
parser.add_argument("--mode_GWAS", action='store_true')
parser.add_argument("--input_path", type=str)
parser.add_argument("--input_path_loci", type=str)
parser.add_argument("--output_path", type=str)
parser.add_argument("--output_path_final", type=str)
#files
parser.add_argument("--hmagma_name_file", type=str)
parser.add_argument("--gene_update_file", type=str)
parser.add_argument("--gene_annotation_file", type=str)
#parameters
parser.add_argument("--consensus_dict", type=str)
parser.add_argument("--coloc_pp4", type=str)


args = parser.parse_args()

mode_pleiotropy = args.mode_pleiotropy
mode_GWAS = args.mode_GWAS
input_path = args.input_path #path of input
input_path_loci = args.input_path_loci #input path for loci parameters
output_path = args.output_path #path of output
output_path_final = args.output_path_final #path for final output
#files
hmagma_name_file = args.hmagma_name_file
gene_update_file = args.gene_update_file
gene_annotation_file = args.gene_annotation_file
#dictionary
consensus_dict = eval(args.consensus_dict)
eqtl_tissue_dict = consensus_dict['eQTL'] #eQTL tissue dictionary
sqtl_tissue_dict = consensus_dict['sQTL'] #sQTL tissue dictionary
pqtl_tissue_dict = consensus_dict['pQTL'] #pQTL tissue dictionary
clinvar_hpo_dict = consensus_dict['ClinVar'] #HPO dictionary
mgi_mpo_dict = consensus_dict['MGI'] #HPO dictionary
hmagma_tissue_dict = consensus_dict['HMAGMA'] #H-MAGMA tissue dictionary
#PP4
coloc_pp4 = eval(args.coloc_pp4)
coloc_eqtl_pp4,coloc_sqtl_pp4,coloc_pqtl_pp4 = coloc_pp4

#paramerte--------------------------------------------

loci = pd.read_csv('%s/loci.parameter' % (input_path_loci),sep=' ',header=None)
loci.columns = ['Prefix','Locus','indsig','chr','s','e']
loci = loci.drop_duplicates(subset=['Prefix','Locus'])

#eQTL
eqtl = pd.read_csv('%s/eQTL.csv' % (input_path))
#Drop invalid rows first
eqtl['PP.H4'] = eqtl['PP.H4'].apply(pd.to_numeric, errors='coerce')
eqtl = eqtl.dropna().drop_duplicates()
#Drop invalid rows first
eqtl = eqtl[eqtl['PP.H4']>=coloc_eqtl_pp4]
eqtl['Gene_id'] = eqtl['Gene_id'].apply(lambda x: x.split('.')[0])

eqtl_result = []
for prefix in loci['Prefix'].unique():
    prefix_tissue = eqtl_tissue_dict.get(prefix,eqtl_tissue_dict.get('all'))
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        eqtl_temp = eqtl[(eqtl['Prefix']==prefix) & (eqtl['Locus']==locus)]
        if len(eqtl_temp) == 0:
            continue
        else:
            for gene in eqtl_temp['Gene'].unique():
                lst = []
                eqtl_temp_gene = eqtl_temp[eqtl_temp['Gene']==gene]
                gene_id = eqtl_temp_gene['Gene_id'].values[0]
                tissue_lst = eqtl_temp_gene['Tissue'].to_list()
                lst = [prefix,locus,gene,gene_id]
                for group in prefix_tissue.keys():
                    lst += [','.join([x for x in tissue_lst if x in prefix_tissue[group].split(',')])]
                eqtl_result.append(lst)
#max group
max_g_eQTL = max([len(eqtl_tissue_dict[prefix]) for prefix in eqtl_tissue_dict.keys()])
eqtl_result_df = pd.DataFrame(eqtl_result,columns=['Prefix','Locus','Gene','Gene_id']+['eQTL_group_%i' % (i) for i in range(1,max_g_eQTL+1)])


#sQTL
sqtl = pd.read_csv('%s/sQTL.csv' % (input_path))
#Drop invalid rows first
sqtl['PP.H4'] = sqtl['PP.H4'].apply(pd.to_numeric, errors='coerce')
sqtl = sqtl.dropna().drop_duplicates()
#Drop invalid rows first
sqtl = sqtl[sqtl['PP.H4']>=coloc_sqtl_pp4]
sqtl['Gene_id'] = sqtl['Gene_id'].apply(lambda x: x.split('.')[0])

sqtl_result = []
for prefix in loci['Prefix'].unique():
    prefix_tissue = sqtl_tissue_dict.get(prefix, sqtl_tissue_dict.get('all'))
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        sqtl_temp = sqtl[(sqtl['Prefix']==prefix) & (sqtl['Locus']==locus)]
        if len(sqtl_temp) == 0:
            continue
        else:
            for gene in sqtl_temp['Gene'].unique():
                lst = []
                sqtl_temp_gene = sqtl_temp[sqtl_temp['Gene']==gene]
                gene_id = sqtl_temp_gene['Gene_id'].values[0]
                tissue_lst = sqtl_temp_gene['Tissue'].to_list() 

                lst = [prefix, locus, gene, gene_id]
                for group in prefix_tissue.keys():
                    filtered_unique_tissues = sorted(list(set([x for x in tissue_lst if x in prefix_tissue[group].split(',')])))
                    lst += [','.join(filtered_unique_tissues)]
                sqtl_result.append(lst)
#max group
max_g_sQTL = max([len(sqtl_tissue_dict[prefix]) for prefix in sqtl_tissue_dict.keys()])
sqtl_result_df = pd.DataFrame(sqtl_result,columns=['Prefix','Locus','Gene','Gene_id']+['sQTL_group_%i' % (i) for i in range(1,max_g_sQTL+1)])


#pQTL
pqtl = pd.read_csv('%s/pQTL.csv' % (input_path))
#Drop invalid rows first
pqtl['PP.H4'] = pqtl['PP.H4'].apply(pd.to_numeric, errors='coerce')
pqtl = pqtl.dropna().drop_duplicates()
#Drop invalid rows first
pqtl = pqtl[pqtl['PP.H4']>=coloc_pqtl_pp4]
pqtl['Gene_id'] = pqtl['Gene_id'].apply(lambda x: x.split('.')[0])

pqtl_result = []
for prefix in loci['Prefix'].unique():
    prefix_tissue = pqtl_tissue_dict.get(prefix,pqtl_tissue_dict.get('all'))
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        pqtl_temp = pqtl[(pqtl['Prefix']==prefix) & (pqtl['Locus']==locus)]
        if len(pqtl_temp) == 0:
            continue
        else:
            for gene in pqtl_temp['Gene'].unique():
                lst = []
                pqtl_temp_gene = pqtl_temp[pqtl_temp['Gene']==gene]
                gene_id = pqtl_temp_gene['Gene_id'].values[0]
                tissue_lst = pqtl_temp_gene['Tissue'].to_list()
                lst = [prefix,locus,gene,gene_id]
                for group in prefix_tissue.keys():
                    lst += [','.join([x for x in tissue_lst if x in prefix_tissue[group].split(',')])]
                pqtl_result.append(lst)
#max group
max_g_pQTL = max([len(pqtl_tissue_dict[prefix]) for prefix in pqtl_tissue_dict.keys()])
pqtl_result_df = pd.DataFrame(pqtl_result,columns=['Prefix','Locus','Gene','Gene_id']+['pQTL_group_%i' % (i) for i in range(1,max_g_pQTL+1)])


#PAV
pav = pd.read_csv('%s/PAV_eligible.csv' % (input_path))
pav_result = []
for prefix in loci['Prefix'].unique():
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        pav_temp = pav[(pav['prefix']==prefix) & (pav['loci']==locus)]
        if len(pav_temp) == 0:
            continue
        else:
            for gene in pav_temp['Gene'].unique():
                pav_temp_gene = pav_temp[pav_temp['Gene']==gene]
                gene_id = pav_temp_gene['Gene_id'].values[0]
                conseq = str(dict(pav_temp_gene[['proxy_snp','Consequence']].values))
                conseq = re.sub("{|}|'", "", conseq)
                pav_result.append([prefix,locus,gene,gene_id,conseq])
pav_result_df = pd.DataFrame(pav_result,columns=['Prefix','Locus','Gene','Gene_id','PAV'])


#ClinVar
clinvar = pd.read_csv('%s/ClinVar.csv' % (input_path))
clinvar_result = []
for prefix in loci['Prefix'].unique():
    loci_prefix = loci[loci['Prefix']==prefix]
    prefix_hpo = clinvar_hpo_dict.get(prefix,clinvar_hpo_dict.get('all'))
    for locus in loci_prefix['Locus'].unique():
        clinvar_temp = clinvar[(clinvar['Prefix']==prefix) & (clinvar['Locus']==locus)]
        if len(clinvar_temp) == 0:
            continue
        else:
            for gene in clinvar_temp['Gene'].unique():
                clinvar_temp_gene = clinvar_temp[clinvar_temp['Gene']==gene]
                gene_id = clinvar_temp_gene['Gene_id'].values[0]
                lst = [prefix,locus,gene,gene_id]
                for group in prefix_hpo.keys():
                    hpo_lst = prefix_hpo[group].split(',')
                    clinvar_temp_gene_ = clinvar_temp_gene[clinvar_temp_gene['HPO_original'].isin(hpo_lst)][['MedGene_name','HPO']].dropna().drop_duplicates()
                    result = clinvar_temp_gene_.groupby('MedGene_name')['HPO'].apply(lambda x: '; '.join(set('; '.join(x).split('; '))))
                    lst += [str({k: v for k, v in result.items()})]
                clinvar_result.append(lst)
max_g_ClinVar = max([len(clinvar_hpo_dict[prefix]) for prefix in clinvar_hpo_dict.keys()])
clinvar_result_df = pd.DataFrame(clinvar_result,columns=['Prefix','Locus','Gene','Gene_id']+['ClinVar_group_%i' % (i) for i in range(1,max_g_ClinVar+1)])


#MGI
mgi = pd.read_csv('%s/MGI.csv' % (input_path))

mgi_result = []
for prefix in loci['Prefix'].unique():
    loci_prefix = loci[loci['Prefix']==prefix]
    prefix_mpo = mgi_mpo_dict.get(prefix,mgi_mpo_dict.get('all'))
    for locus in loci_prefix['Locus'].unique():
        mgi_temp = mgi[(mgi['Prefix']==prefix) & (mgi['Locus']==locus)]
        if len(mgi_temp) == 0:
            continue
        else:
            for gene in mgi_temp['Gene'].unique():
                mgi_temp_gene = mgi_temp[mgi_temp['Gene']==gene]
                gene_id = mgi_temp_gene['Gene_id'].values[0]
                lst = [prefix,locus,gene,gene_id]
                for group in prefix_mpo.keys():
                    mpo_lst = prefix_mpo[group].split(',')
                    mgi_temp_gene_ = mgi_temp_gene[mgi_temp_gene['MPO_original'].isin(mpo_lst)][['Allele_symbol','MPO']].dropna().drop_duplicates()
                    result = mgi_temp_gene_.groupby('Allele_symbol')['MPO'].apply(lambda x: '; '.join(set('; '.join(x).split('; '))))
                    lst += [str({k: v for k, v in result.items()})]
                mgi_result.append(lst)
max_g_MGI = max([len(mgi_mpo_dict[prefix]) for prefix in mgi_mpo_dict.keys()])
mgi_result_df = pd.DataFrame(mgi_result,columns=['Prefix','Locus','Gene','Gene_id']+['MGI_group_%i' % (i) for i in range(1,max_g_MGI+1)])


#H-MAGMA
hmagma = pd.read_csv('%s/H-MAGMA_significant.csv' % (input_path))
hmagma_tissue = pd.read_excel(hmagma_name_file,header=None).dropna()
tissue_dict = {i.replace(' ',''):j.replace(' ','_') for i,j in hmagma_tissue[[0,1]].values}

hmagma_result = []
for prefix in loci['Prefix'].unique():
    prefix_tissue = hmagma_tissue_dict.get(prefix,hmagma_tissue_dict.get('all'))
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        hmagma_temp = hmagma[(hmagma['Prefix']==prefix) & (hmagma['Locus']==locus)]
        if len(hmagma_temp) == 0:
            continue
        else:
            for gene in hmagma_temp['Gene'].unique():
                lst = []
                hmagma_temp_gene = hmagma_temp[hmagma_temp['Gene']==gene]
                gene_id = hmagma_temp_gene['Gene_id'].values[0]
                tissue_lst = hmagma_temp_gene['Tissue'].to_list()
                lst = [prefix,locus,gene,gene_id]
                for group in prefix_tissue.keys():
                    lst += [','.join([tissue_dict[x] for x in tissue_lst if x in prefix_tissue[group].split(',')])]
                hmagma_result.append(lst)
#max group
max_g_HMAGMA = max([len(hmagma_tissue_dict[prefix]) for prefix in hmagma_tissue_dict.keys()])
hmagma_result_df = pd.DataFrame(hmagma_result,columns=['Prefix','Locus','Gene','Gene_id']+['HMAGMA_group_%i' % (i) for i in range(1,max_g_HMAGMA+1)])


#PoPS
pops = pd.read_csv('%s/PoPS.csv' % (input_path))
pops_result = []
for prefix in loci['Prefix'].unique():
    loci_prefix = loci[loci['Prefix']==prefix]
    for locus in loci_prefix['Locus'].unique():
        pops_temp = pops[(pops['Prefix']==prefix) & (pops['Locus']==locus)]
        if len(pops_temp) == 0:
            continue
        else:
            for gene in pops_temp['Gene'].unique():
                pops_temp_gene = pops_temp[pops_temp['Gene']==gene]
                gene_id = pops_temp_gene['Gene_id'].values[0]
                rank = pops_temp_gene['Rank'].values[0]
                pops_result.append([prefix,locus,gene,gene_id,rank])
pops_result_df = pd.DataFrame(pops_result,columns=['Prefix','Locus','Gene','Gene_id','pops_rank'])



#merge the evidence
gene_update = pd.read_csv(gene_update_file,sep='\t')
gene_update_dict = dict(gene_update[['Requested ID','Matched ID(s)']].values)
gene_name = pd.read_csv(gene_annotation_file,sep='\t')
id_name_dict = dict(gene_name[['Gene_ID','NAME']].values)

#update ID
for df in [eqtl_result_df,sqtl_result_df,pqtl_result_df,pav_result_df,clinvar_result_df,mgi_result_df,hmagma_result_df,pops_result_df]:
    df['Gene_id'] = df['Gene_id'].apply(lambda x: gene_update_dict.get(x,x))
    df['Gene'] = df.apply(lambda row: id_name_dict.get(row['Gene_id'],row['Gene']),axis=1)

#merge
merged = []
for df in [eqtl_result_df,sqtl_result_df,pqtl_result_df,pav_result_df,clinvar_result_df,mgi_result_df,hmagma_result_df,pops_result_df]:
    try:
        merged = pd.merge(merged,df,on=['Prefix','Locus','Gene_id','Gene'],how='outer')
    except:
        merged = df.copy()
merged = merged.mask(merged == '') #replace the '' value
merged = merged.mask(merged == "{}") #replace the {} value

#score for each group
for module in ['eQTL','sQTL','pQTL','ClinVar','MGI','HMAGMA']:
    n_group = eval('max_g_%s' % module)
    for i in range(1,n_group+1):
        merged['%s_score_g%i' % (module,i)] = merged['%s_group_%i' % (module,i)].apply(lambda x: 0 if pd.isna(x) else 1)
#PAV
merged['PAV_score'] = 1 - merged['PAV'].apply(lambda x: pd.isna(x))
#PoPS
merged['PoPS_score'] = merged['pops_rank'].apply(lambda x: 1 if x<=3 else 0)

if mode_pleiotropy:
    for module in ['eQTL','sQTL','pQTL','ClinVar','MGI','HMAGMA']:
        n_group = eval('max_g_%s' % module)
        if n_group==2:
            merged['%s_score' % (module)] = merged[['%s_score_g1' % (module),'%s_score_g2' % (module)]].min(axis=1)
        elif n_group>2:
            extra_columns_n = n_group - 2
            extra_columns = ['%s_score_g%i' % (module,i) for i in range(3,extra_columns_n)]
            merged['%s_score' % (module)] = merged[['%s_score_g1' % (module),'%s_score_g2' % (module)]].min(axis=1) + merged[extra_columns].sum(axis=1)
        elif n_group==1 and module == 'pQTL':
            merged['%s_score' % (module)] = merged['%s_score_g1' % (module)]
        else:
            raise ValueError('Number of groups <2 for %s' % (module))
if mode_GWAS:    
    for module in ['eQTL','sQTL','pQTL','ClinVar','MGI','HMAGMA']:
        n_group = eval('max_g_%s' % module)
        columns = ['%s_score_g%i' % (module,i) for i in range(1,n_group+1)]
        merged['%s_score' % (module)] = merged[columns].sum(axis=1)
       
#total score
merged['Total_score'] = merged[['eQTL_score','sQTL_score','pQTL_score','PAV_score','ClinVar_score','MGI_score',
                                'HMAGMA_score','PoPS_score']].sum(axis=1)
merged = merged[merged['Total_score']>=1]

eqtl_result_df.to_csv('%s/eQTL.csv' % (output_path),index=False)
sqtl_result_df.to_csv('%s/sQTL.csv' % (output_path),index=False)
pqtl_result_df.to_csv('%s/pQTL.csv' % (output_path),index=False)
pav_result_df.to_csv('%s/PAV.csv' % (output_path),index=False)
clinvar_result_df.to_csv('%s/ClinVar.csv' % (output_path),index=False)
mgi_result_df.to_csv('%s/MGI.csv' % (output_path),index=False)
hmagma_result_df.to_csv('%s/H-MAGMA.csv' % (output_path),index=False)
pops_result_df.to_csv('%s/PoPS.csv' % (output_path),index=False)
merged = merged.sort_values(by=['Prefix','Locus','Total_score'],ascending=False)
merged.to_csv('%s/Consensus.csv' % (output_path_final),index=False)

mapped = merged[merged['Total_score']>=3].copy()
if len(mapped) > 0:
    mapped = mapped[mapped['Total_score'] == mapped.groupby(['Prefix','Locus'])['Total_score'].transform('max')]
mapped.to_csv('%s/Consensus_mapped_genes.csv' % (output_path_final),index=False)





















