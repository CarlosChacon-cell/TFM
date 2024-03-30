
#This is just a small code to compare the pae interaction of the hits with the values obtained through the pymol detected interactions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pymol_scores=pd.read_csv('/emdata/cchacon/hits_compilation_20240227_NGR/all/scores.csv', header = 0, sep ='\t')
print(pymol_scores)

af2_scores=pd.read_csv('/emdata/cchacon/hits_compilation_20240227_NGR/all.csv', header=0, sep=',')
print(af2_scores)

df=pd.merge(pymol_scores, af2_scores, on='description')
print(df)

pattern=r'^run_\d+'
for i in range(0,len(df['description'])):
    descriptor=df['description'][i]
    target_name= re.search(pattern, descriptor)
    target_name=target_name.group(0)
    df['description'][i]=target_name



ax=df.plot.scatter(x='SCORE', y='pae_interaction')

for idx, row in df.iterrows():
    ax.annotate(row['description'], (row['SCORE'], row['pae_interaction']))
