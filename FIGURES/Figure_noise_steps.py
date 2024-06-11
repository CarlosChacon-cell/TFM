import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import matplotlib.patches as patches
import scipy
from statannotations.Annotator import Annotator
from scipy.stats import ks_2samp
# Load data
df = pd.read_csv('/emdata/cchacon/RFD_PD_test/pdbs/rmsd_noise_steps.csv')

# Prepare data for plotting
df_filtered = df[df['campaign'].str.contains('campaign_run_112_ns')]
df_filtered['Noise'] = df_filtered['campaign'].str.extract('campaign_run_112_ns(\d+)').astype(int)

df_whole=pd.read_csv('/emdata/cchacon/RFD_PD_test/output_af2_NS.csv')
print(df_whole)

df_barplot=pd.DataFrame()
df_barplot['Noise']=df_whole['Noise_Steps'].unique()
percentage=[]
# Loop through each unique noise value
for noise in df_whole['Noise_Steps'].unique():
    # Filter the DataFrame based on the current noise level and the conditions
    filtered = df_whole[(df_whole['Noise_Steps'] == noise) & 
                        (df_whole['pae_interaction'] < 10) & 
                        (df_whole['plddt_binder'] > 80)]
    # Calculate the percentage
    percent = len(filtered) / len(df_whole[df_whole['Noise_Steps'] == noise])
    
    # Append the result to the percentage list
    percentage.append(percent)
print(percentage)
df_barplot['Percentage']=percentage
#Barplot success rate
plt.figure(figsize=(10, 12))
barplot = sns.barplot(x='Noise', y='Percentage', data=df_barplot, palette='crest')

# Add titles and labels
plt.title('Noise Steps Success Rates', fontsize=20)
plt.ylabel('Hits Success Rate (%)', fontsize=16)

plt.xlabel('Noise Steps (#)', fontsize=16)
barplot.set_xticklabels((f'5\n n={len( df_whole[df_whole["Noise_Steps"]==5])}', 
                         f'10\n n={len(df_whole[df_whole["Noise_Steps"]==10])}', 
                         f'20\n n={len(df_whole[df_whole["Noise_Steps"]==20])}', 
                         f'30\n n={len(df_whole[df_whole["Noise_Steps"]==30])}'), fontsize=14)
plt.tick_params(axis='y', labelsize=14)
# Show the plot
plt.ylim((0,0.07))
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/barplot_noisesteps.png')








X=[5,10,20,30]
y=[]
for noise in df_filtered['Noise'].unique():
    y.append(np.mean(df_filtered['RMSD'][df_filtered['Noise']==noise]))

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
y_pred = [intercept + slope * value for value in X]

# Initialize the plot

plt.figure(figsize=(10,12))
sns.set_theme(style='whitegrid')
sns.pointplot(data=df_filtered, x='Noise', y='RMSD', hue= 'Noise', palette='crest', s=100,capsize=0.2, ci=95)


# Plot the regression line
# sns.lineplot(x=X, y=y_pred, label=f'$r$={r_value:.2f}', color='red', ax=ax) ##Doesn't work##
plt.title('RMSD vs Noise Steps', fontsize=20)
plt.xlabel('Noise Steps (#)', fontsize=16)
plt.ylabel(f'RMSD ($\AA$)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)
plt.legend(title='Noise Steps')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseStepsRMSD.png')


df_filtered=df_filtered[df_filtered['Noise']!=30]
plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_filtered, x='Noise', y='RMSD', palette='crest')
annotator=Annotator(
    violinplot,
    data=df_filtered,
    x='Noise',
    y='RMSD',
    pairs=[
        (5,10),(5,20),
        (10,20)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('RMSD vs Noise Steps', fontsize=20)
plt.ylabel(f'RMSD ($\AA$)', fontsize=16)
plt.xlabel('Noise Steps(#)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_rmsd_noisesteps.png')

#df_filter hits
df_hits=df_filtered[(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder']>80)]
d_stat,p_value=ks_2samp(df_hits['RMSD'][df_hits['Noise']==10], df_hits['RMSD'][df_hits['Noise']==5])
print('P_VALUE BY KS IS: ',  p_value)

#A violinplot in case it is interesting 
plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_hits, x='Noise', y='pae_interaction', palette='crest')
annotator=Annotator(
    violinplot,
    data=df_hits,
    x='Noise',
    y='pae_interaction',
    pairs=[
        (5,10),(5,20),
        (10,20)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Steps vs PAE interaction', fontsize=20)
plt.xlabel('Noise Steps (#)', fontsize=16)
plt.ylabel(f'PAE interaction ($\AA$)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)


plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_pae_interaction_noise_steps.png')

plt.figure(figsize=(10,12))

violinplot=sns.violinplot(data=df_hits, x='Noise', y='plddt_binder', palette='crest')
violinplot.set_ylim((78,88))

annotator=Annotator(
    violinplot,
    data=df_hits,
    x='Noise',
    y='plddt_binder',
    pairs=[
        (5,10),(5,20),
        (10,20)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Steps vs pLDDT binder', fontsize=20)
plt.xlabel('Noise Steps (#)', fontsize=16)
plt.ylabel(f'pLDDT binder (Arb.Units)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)


plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_plddt_binder_noise_steps.png')



# Load the dataset
df = pd.read_csv('/emdata/cchacon/20240522_run_264_noise_steps/output_af2.csv')

df_filtered=df[df['campaign'].str.contains('campaign')]
df_filtered['Noise']=df_filtered['campaign'].str.extract('campaign_(\d+)ns').astype(int)
df_5ns=df_filtered[df_filtered['Noise']==5]
df_10ns=df_filtered[df_filtered['Noise']==10]
df_20ns=df_filtered[df_filtered['Noise']==20]
df_30ns=df_filtered[df_filtered['Noise']==30]

# Original data point
original = {
    'pae_interaction': [7.383],
    'plddt_binder': [88.05]
}

# Convert the original data to a DataFrame
df_original = pd.DataFrame(original)

# Set the style
sns.set(style="whitegrid")

# Create a scatter plot
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df_filtered, x='plddt_binder', y='pae_interaction', hue='Noise', s=100, palette='crest')
sns.scatterplot(data=df_original, x='plddt_binder', y='pae_interaction', color='red', s=100, label='Original')

# Add titles and labels
plt.title('PAE Interaction vs. pLDDT Binder run_264', fontsize=20)
plt.xlabel('PAE Interaction', fontsize=16)
plt.ylabel('pLDDT Binder', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)


# Highlight specific region with a square
rect = patches.Rectangle((0, 80), 10, 20, linewidth=2, edgecolor='black', facecolor='none')
plt.gca().add_patch(rect)
plt.text(4, 81, 'Hits Zone', color='black', fontsize=12)
# Invert the x-axis
plt.gca().invert_xaxis()
# Add legend
plt.legend(title='Noise Steps')
# Add grid
plt.grid(True)

# Show the plot
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseStepsRun264Figure.png')

