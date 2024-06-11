
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy

all_df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv', header=0)
filtered_df=all_df[(all_df['pae_interaction'] < 10) & (all_df['plddt_binder'] > 80)]
filtered_df=filtered_df[filtered_df['campaign'].str.contains('run_112_noise_scale')]
filtered_df['Noise'] = filtered_df['campaign'].str.extract('run_112_noise_scale_(\d+)').astype(str)
filtered_df.loc[filtered_df['Noise'] == '005', 'Noise'] = 0.05
filtered_df.loc[filtered_df['Noise'] == '01', 'Noise'] = 0.1
filtered_df.loc[filtered_df['Noise'] == '05', 'Noise'] = 0.5
filtered_df.loc[filtered_df['Noise'] == '1', 'Noise'] = 1.0

nw_df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/pdbs/pairwise_alignment.csv', header=0, names=['score','identities','description','target'])
merge_df=pd.merge(nw_df,filtered_df, on='description', how='inner')

x=merge_df['Noise'].unique().tolist()
y=[]
for noise in merge_df['Noise'].unique():
    y.append(np.mean(merge_df['score'][merge_df['Noise']==noise]))
print(x)
print(y)

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

print('R_VALUE IS: ', r_value)


plt.figure(figsize=(12,8))
sns.pointplot(data=merge_df, x='Noise', y='score', hue='Noise', palette='crest', s=100, capsize=0.2, ci=95)
plt.title('NW score vs Noise Scale')
plt.xlabel('Noise Scale (Arb.U)')
plt.ylabel('Score')
plt.legend(title='Noise Scale')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NWscore_noisescale.png')
