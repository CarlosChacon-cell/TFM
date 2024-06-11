
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy



all_df=pd.read_csv('/emdata/cchacon/RFD_PD_test/output_af2.csv', header=0)
filtered_df=all_df[(all_df['pae_interaction'] < 10) & (all_df['plddt_binder'] > 80)]
filtered_df=filtered_df[filtered_df['campaign'].str.contains('campaign_run_112_ns')]
filtered_df['Noise'] = filtered_df['campaign'].str.extract('campaign_run_112_ns(\d+)').astype(int)

nw_df=pd.read_csv('/emdata/cchacon/RFD_PD_test/pdbs/pairwise_alignment.csv', header=0, names=['score','identities','description','target'])
merge_df=pd.merge(nw_df,filtered_df, on='description', how='inner')

x=merge_df['Noise'].unique().tolist()
y=[]
for noise in merge_df['Noise'].unique():
    y.append(np.mean(merge_df['score'][merge_df['Noise']==noise]))


slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

print('R_VALUE IS: ', r_value)


plt.figure(figsize=(12,8))
sns.pointplot(data=merge_df, x='Noise', y='score', hue='Noise', palette='crest', s=100, capsize=0.2, ci=95)
plt.title('NW score vs Noise')
plt.xlabel('Noise Steps (#)')
plt.legend(title='Noise Steps')
plt.ylabel('Score')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NWscore_noisesteps.png')
