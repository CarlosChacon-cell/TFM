
#This is just a small code to compare the pae interaction of the hits with the values obtained through the pymol detected interactions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re 

pymol_scores=pd.read_csv('/emdata/cchacon/hits_compilation_20240227_NGR/all/scores.csv', header = 0, sep =',')
print(pymol_scores.keys())

af2_scores=pd.read_csv('/emdata/cchacon/hits_compilation_20240227_NGR/all.csv', header=0, sep=',')
stripped_keys = {key.strip(): value for key, value in af2_scores.items()}
af2_scores = pd.DataFrame(stripped_keys)
print(af2_scores.keys())

df=pd.merge(pymol_scores, af2_scores, on='description')
df.to_csv('plot.csv', index=False)



print(df)

pattern=r'^run_\d+'
for i in range(0,len(df['description'])):
    descriptor=df['description'][i]
    target_name= re.search(pattern, descriptor)
    try:
        target_name=target_name.group(0)
    except AttributeError:
        pattern1=r'alpha_\d+'
        target_name= re.search(pattern1,descriptor)
        target_name=target_name.group(0)
    df['description'][i]=target_name



ax=df.plot.scatter(x='SCORE', y='pae_interaction')

for idx, row in df.iterrows():
    ax.annotate(row['description'], (row['SCORE'], row['pae_interaction']))
