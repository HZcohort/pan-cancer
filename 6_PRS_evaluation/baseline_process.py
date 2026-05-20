from pathlib import Path

import pandas as pd
import numpy as np

def drinking_cate (status,sex,eid,dict_):
    if status < 0 or pd.isna(status):
        return '6-nan'
    elif status == 1:
        return '4-previous'
    elif status == 0:
        return '5-never'
    elif status == 2:
        gday = dict_[eid]
        if pd.isna(gday):
            return '6-nan'
        elif gday <= 14:
            return '1-within'
        elif gday <= 35 and sex == 0:
            return '2-hazardous'
        elif gday <= 50 and sex == 1:
            return '2-hazardous'
        elif gday > 35 and sex == 0:
            return '3-harmful'
        elif gday > 50 and sex == 1:
            return '3-harmful'
        else:
            return 'error'
    else:
        return 'error'

def smoking_cate (status,eid,dict_):
    if status < 0 or pd.isna(status):
        return '6-nan'
    elif status == 1:
        return '4-previous'
    elif status == 0:
        return '5-never'
    elif status == 2:
        numday = dict_[eid]
        if pd.isna(numday):
            return '6-nan'
        elif numday < 5:
            return '1-within'
        elif numday < 15:
            return '2-hazardous'
        elif numday >= 15:
            return '3-harmful'
        else:
            return 'error'
    else:
        return 'error'

BASE_DIR = Path("/path/to/ukb_prs_evaluation")
input_dir = BASE_DIR / "input"
BASE_DIR.mkdir(parents=True, exist_ok=True)
output_path = BASE_DIR / "baseline.csv"

baseline = pd.read_csv(input_dir / "baseline.csv",index_col=0)
new_name = {'sex_f31_0_0':'sex', 'year_of_birth_f34_0_0':'birth_year', 'month_of_birth_f52_0_0':'birth_month',
            'townsend_deprivation_index_at_recruitment_f189_0_0':'social','body_mass_index_bmi_f21001_0_0':'BMI',
            'date_of_attending_assessment_centre_f53_0_0':'date_attend', 'tea_intake_f1488_0_0':'tea_intake',
            'coffee_type_f1508_0_0':'coffee_type', 'smoking_status_f20116_0_0':'smoking',
            'alcohol_drinker_status_f20117_0_0':'drinking', 'coffee_intake_f1498_0_0':'coffee_intake',
            'average_total_household_income_before_tax_f738_0_0':'income',
            'summed_met_minutes_per_week_for_all_activity_f22040_0_0':'physical',
            'cooked_vegetable_intake_f1289_0_0':'vegetable-1','salad_raw_vegetable_intake_f1299_0_0':'vegetable-2',
            'fresh_fruit_intake_f1309_0_0':'fruit-1','dried_fruit_intake_f1319_0_0':'fruit-2'}
baseline = baseline.rename(columns=new_name)
baseline = baseline.loc[baseline['eid']!=6025218]
len_ori = len(baseline)

withdraw = pd.read_csv(input_dir / "ukb_eid.csv.gz")
baseline = baseline[baseline['eid'].isin(withdraw['Participant ID'])]

baseline['social_cate'] = pd.cut(baseline['social'],[-np.inf,-3.64,-2.14,0.55,np.inf],right=False,
                                 labels=['1-<-3.64','2--3.64--2.14','3--2.14-0.55','4->=0.55']).values.add_categories('5-nan')
baseline['social_cate'] = baseline['social_cate'].fillna('5-nan')
baseline['birth_date'] = baseline.apply(lambda row: f"{int(row['birth_year'])}-{int(row['birth_month'])}-01",axis=1)

baseline['age'] = (pd.to_datetime(baseline['date_attend']) - pd.to_datetime(baseline['birth_date'])).dt.days/365.25
baseline['BMI'] = pd.cut(baseline['BMI'],[-np.inf,24.1,29.9,np.inf],right=False,
                         labels=['1-<24.1','2-24.1-29.9','3->=29.9']).values.add_categories('4-nan')
baseline['BMI'] = baseline['BMI'].fillna('4-nan')

#merge with alcohol drinking
g_day = pd.read_csv(input_dir / "lifestyle.csv",usecols=['eid','alcohol_gday'])
g_day['alcohol_gday'] = g_day['alcohol_gday'].apply(lambda x: x/8*7)
g_day_dict = {int(i):j for i,j in g_day.values}
baseline['drinking'] = baseline.apply(lambda row: drinking_cate(row['drinking'],row['sex'],row['eid'],g_day_dict),axis=1)
#merge with alcohol drinking

#merge with smoking amount
smoking_num = pd.read_csv(input_dir / "lifestyle.csv",usecols=['smoking_num','eid'])
smoking_num_dict = {int(i):j for i,j in smoking_num.values}
baseline['smoking'] = baseline.apply(lambda row: smoking_cate(row['smoking'],row['eid'],smoking_num_dict),axis=1)
#merge with smoking amount

baseline['income'] = baseline['income'].apply(lambda x: '1-<18000' if x==1
                                                    else '2-18000-51999' if x==2 or x==3
                                                    else '3->52000' if x>=4
                                                    else '4-nan')

baseline['physical'] = pd.cut(baseline['physical'],[-np.inf,798,3552,np.inf],right=False,
                              labels=['1-<798','2-798-3552','3->=3552']).values.add_categories('4-nan')
baseline['physical'] = baseline['physical'].fillna('4-nan')

##diet
for var in ['fruit-1','fruit-2','vegetable-1','vegetable-2']:
    baseline[var] = baseline[var].apply(lambda x: 0.5 if x==-10 else x if x>=0 else np.nan)
diet = baseline.dropna(subset=['fruit-1','fruit-2','vegetable-1','vegetable-2'],how='all')
diet['diet'] = np.nansum(diet[['fruit-1','fruit-2']].values,axis=1) + np.nansum(diet[['vegetable-1','vegetable-2']].values,axis=1)*21.25/80
baseline = pd.merge(baseline,diet[['eid','diet']],on='eid',how='left')
baseline['diet'] = baseline['diet'].apply(lambda x: '1-notenough' if x<5 else '2-enough' if x>=5 else '3-nan')

#merge with PCA data
geno_qc = pd.read_csv(input_dir / "genotype_principal_components.csv",
                      usecols=['eid','genetic_principal_components_f22009_0_1',
                               'genetic_principal_components_f22009_0_2','genetic_principal_components_f22009_0_3',
                               'genetic_principal_components_f22009_0_4','genetic_principal_components_f22009_0_5'])
geno_qc = geno_qc.rename(columns={'genetic_principal_components_f22009_0_1':'pca1','genetic_principal_components_f22009_0_2':'pca2',
                                  'genetic_principal_components_f22009_0_3':'pca3','genetic_principal_components_f22009_0_4':'pca4',
                                  'genetic_principal_components_f22009_0_5':'pca5'}).dropna()
baseline = pd.merge(baseline,geno_qc,on='eid',how='inner')

#sex code
baseline['sex'] = baseline['sex'].map({1:0,0:1})

#date of death/censor
censor = pd.read_csv(input_dir / "data_death_censor.csv")
censor = censor.rename(columns={'Participant ID':'eid','Date of censor':'censor_date'})
censor['censor_date'] = pd.to_datetime(censor['censor_date'])
censor = censor.sort_values(by='censor_date')
censor = censor.drop_duplicates(subset=['eid'],keep='first')
baseline = pd.merge(baseline,censor[['eid','censor_date']],how='left',on='eid')
baseline['date_end'] = baseline['censor_date'].fillna(pd.to_datetime('2022-10-31'))

baseline.to_csv(output_path,index=False)

