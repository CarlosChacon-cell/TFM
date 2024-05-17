
import pandas as pd 

df=pd.read_csv('/emdata/cchacon/inmuno/micropeptides/run_33_europe.txt', sep='\t', skiprows=11)
df_subset = df.iloc[:, :5] 
df_filtered=df_subset[(df_subset['%Rank_bestAllele'] < 0.5) & (df_subset['%RankBinding_bestAllele'] < 0.5)]
print(df_filtered)


