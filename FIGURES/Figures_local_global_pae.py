import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns 
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import ttest_ind
from statannotations.Annotator import Annotator
from sklearn.metrics import r2_score
# Load the CSV file
file_path = glob.glob('/emdata/cchacon/Scaffolding/campaign_run_9_again/hits/new_cutre/pdbs/pae_local_global.csv')[0]

df = pd.read_csv(file_path)
df_af2=pd.read_csv('/emdata/cchacon/Scaffolding/campaign_run_9_again/hits.csv', header=0)
df_plddt=df_af2[[ 'plddt_binder','description']]
df_plddt.columns=['plddt_binder', 'protein_name']
# df=pd.merge(df,df_plddt, on='protein_name')
# Extract global and local PAE interactions
global_pae = df['pae_interaction_global']
local_pae = df['pae_interaction_local']

# Calculate the differences and improvements
counter = 0
total = 0
differences = []
for i in range(len(df['protein_name'])):
    if df['pae_interaction_local'][i] >= df['pae_interaction_global'][i]:
        differences.append(df['pae_interaction_local'][i] - df['pae_interaction_global'][i])
        counter += 1
        total += 1
    else:
        differences.append(df['pae_interaction_local'][i] - df['pae_interaction_global'][i])
        total += 1

improve = 100 - (counter / total) * 100 

print('Local pae_interaction predicts an improvement for: ', improve, ' of the designs')
print('')
print('Local pae_interaction - global_pae_interaction: ', np.mean(differences))
print('Mean interacting surface: ', np.mean(df['interacting_surface']))

# Initialize the plot
plt.figure(figsize=(10, 14))

# Plot with seaborn
sns.scatterplot(data=df, y='pae_interaction_local', x='pae_interaction_global', hue='interacting_surface', palette='crest', s=100)
plt.xlabel(f'PAE Interaction ($\AA$)', fontsize=22)
plt.ylabel(f'CUTRE ($\AA$)', fontsize=22)
plt.title('Global PAE Interaction vs CUTRE', fontsize=26)
plt.grid(True)
plt.xlim((0, 32))
plt.ylim((0, 32))
# sns.regplot(data=df, y='pae_interaction_local', x='pae_interaction_global', scatter=False, color='black', line_kws={'label': 'Linear fit'})

#R square
slope, intercept, r_value, p_value, std_err = linregress(global_pae, local_pae)
r_squared = r_value**2

print(f'R-squared: {r_squared:.2f}')

# Add vertical and horizontal lines
plt.axvline(x=10, color='red', linestyle='--', label='Global threshold')
plt.axhline(y=10, color='blue', linestyle='--', label='CUTRE threshold')
plt.legend(title='Interacting Surface (%)', title_fontsize=20, fontsize=18 )
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/Figure_newCUTREvsPaeInteraction.png')

# # Define the exponential function
# def exponential_model(x, a, b, c):
#     return a * np.exp(b * x) + c

# # Fit the exponential model
# params_exp, covariance_exp = curve_fit(exponential_model, global_pae, local_pae, maxfev=5000)

# # Generate predictions
# pred_exp = exponential_model(global_pae, *params_exp)

# r2_exp = r2_score(local_pae, pred_exp)
# rss_exp = np.sum((local_pae - pred_exp) ** 2)

# print(f'EXPONENTIAL FITTING R_VALUE: {r2_exp} RSS: {rss_exp}')


# plt.figure(figsize=(10,12))
# no_hit_df=df[(df['pae_interaction_global'] > 20)]
# hit_df=df[(df['pae_interaction_global'] < 20)]
# print(no_hit_df)
# sns.scatterplot(data=no_hit_df, x='pae_interaction_global', y='plddt_binder', color='green', s=100, alpha=0.3)
# sns.scatterplot(data=hit_df, x='pae_interaction_global', y='plddt_binder', hue='protein_name', s=100)
# sns.scatterplot(data=no_hit_df, x='pae_interaction_local', y='plddt_binder', color='blue', s=100, alpha=0.3, marker='s')
# sns.scatterplot(data=hit_df, x='pae_interaction_local', y='plddt_binder', hue='protein_name', s=100, marker='s')
# plt.xlabel('pae_interaction')
# plt.ylabel('plddt_binder')
# plt.gca().invert_xaxis()
# plt.axvline(x=10, color='k', linestyle='-')
# plt.axhline(y=80, color='k', linestyle='-')
# plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/scatterplot_local_vs_global_new_cutre.png')




plt.figure(figsize=(10,14))


#Create the new variables for plotting
df['Improvement']=df['pae_interaction_global']-df['pae_interaction_local']
df['interaction_group']=df['interacting_surface'].round(-1)
hits_0= len(df[df['interaction_group']==0.0])
hits_10=len(df[df['interaction_group']==10.0])
hits_20=len(df[df['interaction_group']==20.0])
hits_30=len(df[df['interaction_group']==30.0])
hits_40=len(df[df['interaction_group']==40.0])


violinplot=sns.boxplot(data=df, x='interaction_group', y='Improvement', palette='Set2')
violinplot.set_xticklabels((f'0\n n={hits_0}', f'10\n n={hits_10}', f'20\n n={hits_20}',f'30\n n={hits_30}',f'40\n n={hits_40}'), fontsize=20)
annotator=Annotator(
    violinplot,
    data=df,
    x='interaction_group',
    y='Improvement',
    pairs=[
        (0,10),(0,20),(0,30),(0,40),
        (10,20),(10,30),(10,40),
        (20,30),(20,40),
        (30,40)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="holm-bonferroni",
    fontsize=18)
    
annotator.apply_and_annotate()
plt.xlabel('Interacting Surface (%)', fontsize=22)
plt.ylabel(f'Difference ($\AA$)', fontsize=22)
plt.title('Metrics Difference vs Interacting Surface', fontsize=26)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)
# sns.scatterplot(data=df, x='interacting_surface', y='Improvement', hue='pae_interaction_global', s=100)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/boxplot_interactingSurface_Improvement_newcutre.png')

# subset_0=df['Improvement'][df['interaction_group']==0]
# subset_10=df['Improvement'][df['interaction_group']==10]
# subset_20=df['Improvement'][df['interaction_group']==20]
# subset_30=df['Improvement'][df['interaction_group']==30]
# subset_40=df['Improvement'][df['interaction_group']==40]

# holm-bonferroni_correction=0.05/10
# t_stat, p_value = ttest_ind(subset_0, subset_10, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('The difference between 0 and 10 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_0, subset_20, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 0 and 20 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_0, subset_30, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference is 0 and 30 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_0, subset_40, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference 0 and 40 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_10, subset_20, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 10 and 20 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_10, subset_30, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 10 amd 30 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_10, subset_40, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 10 amd 40 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_20, subset_30, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 20 amd 30 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_20, subset_40, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 20 amd 40 is significant\n')
#     print('#########\n')
# t_stat, p_value = ttest_ind(subset_40, subset_30, equal_var=False)
# if p_value < holm-bonferroni_correction:
#     print('########\n')
#     print('This difference between 40 and 30 is significant\n')
#     print('#########\n')



