
import pandas as pd 
import glob 


files=glob.glob('./outputs/*txt')
dict_alleles={
    'file':[],
    'Visibility':[]
}
for file in files:
    print(file)
    try:
        df=pd.read_csv(file, sep='\t', skiprows=11)
    except:
        continue
    df_subset = df.iloc[:, :5] 
    df_filtered=df_subset[(df_subset['%Rank_bestAllele'] < 0.5) & (df_subset['%RankBinding_bestAllele'] < 0.5)]
    dict_alleles['file'].append(file.split('/')[2].split('.')[0])
    dict_alleles['Visibility'].append(len(df_filtered))

df_alleles=pd.DataFrame(dict_alleles)
df_alleles.to_csv('Visibility.csv')


