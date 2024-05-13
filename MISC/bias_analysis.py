
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np


df_unbiased=pd.read_csv('antigenic_values_unbiased.csv', index_col=None)
df_wrong=pd.read_csv('antigenic_values_wrong_bias.csv', index_col=None)
df_bias=pd.read_csv('antigenic_values_biased_good.csv', index_col=None)

merged_df=pd.merge(df_unbiased, df_bias, on='protein_name')
merged_df=pd.merge(df_wrong,merged_df, on='protein_name')
merged_df=merged_df.rename(columns={
    'mean_antigenic_x':'Unbiased',
    'mean_antigenic':'+Antigenic',
    'mean_antigenic_y':'-Antigenic'
})



plt.figure(figsize=(8, 6))
# Scatter plot for local PAE with colors based on protein names
plt.scatter(np.zeros(len(merged_df['protein_name'])), merged_df['+Antigenic'], label='+Antigenic',s=100)  # Increased point size (s=100)

# Scatter plot for global PAE with colors based on protein names
plt.scatter(np.ones(len(merged_df['protein_name'])), merged_df['Unbiased'], label='Unbiased', s=100)  # Increased point size (s=100)

# Scatter plot for local PAE with colors based on protein names
plt.scatter(np.full((len(merged_df['protein_name'])),2), merged_df['-Antigenic'], label='-Antigenic',s=100)  # Increased point size (s=100)

for i in range(len(merged_df)):
    plt.plot([0, 1], [merged_df['+Antigenic'], merged_df['Unbiased']], c='k', lw=1)  # Thinner line (lw=1)
    plt.plot([1, 2], [merged_df['Unbiased'], merged_df['-Antigenic']], c='k', lw=1)  # Thinner line (lw=1)

plt.xticks([0, 1, 2], ['+Antigenic', 'Unbiased', '-Antigenic'])
plt.ylabel('Antigenicity')
plt.title('Antigenicity with different bias')
plt.grid(True)
plt.show()

antigenicity_increase=np.mean(merged_df['+Antigenic']-merged_df['Unbiased'])
antigenicity_decrease=np.mean(merged_df['-Antigenic']-merged_df['Unbiased'])

print('Using a antigenic bias has increased the value of antigenicity a mean of: ', antigenicity_increase, '\n')
print('Using an anti-antigenic bias (a genic bias I guess) has decreased the value of antigenicity a mean of: ', antigenicity_decrease, '\n')