import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import scipy
import matplotlib.patches as patches
from statannotations.Annotator import Annotator

# Load data
df = pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv')
print(df)
# Prepare data for plotting
df_filtered_112 = df[df['campaign'].str.contains('run_112_noise_scale')]
df_filtered_112['Noise'] = df_filtered_112['campaign'].str.extract('run_112_noise_scale_(\d+)').astype(str)
df_filtered_112.loc[df_filtered_112['Noise'] == '005', 'Noise'] = 0.05
df_filtered_112.loc[df_filtered_112['Noise'] == '01', 'Noise'] = 0.1
df_filtered_112.loc[df_filtered_112['Noise'] == '05', 'Noise'] = 0.5
df_filtered_112.loc[df_filtered_112['Noise'] == '1', 'Noise'] = 1.0

df_barplot=pd.DataFrame()
df_barplot['Noise']=df_filtered_112['Noise'].unique()
percentage=[]
# Loop through each unique noise value
for noise in df_filtered_112['Noise'].unique():
    # Filter the DataFrame based on the current noise level and the conditions
    filtered = df_filtered_112[(df_filtered_112['Noise'] == noise) & 
                        (df_filtered_112['pae_interaction'] < 10) & 
                        (df_filtered_112['plddt_binder'] > 80)]
    # Calculate the percentage
    percent = len(filtered) / len(df_filtered_112[df_filtered_112['Noise'] == noise])
    
    # Append the result to the percentage list
    percentage.append(percent)
df_barplot['Percentage']=percentage
#Barplot success rate
plt.figure(figsize=(10, 8))

barplot = sns.barplot(x='Noise', y='Percentage', data=df_barplot, palette='crest')

# Add titles and labels
plt.title('Noise Scale Success Rates',fontsize=20)
plt.ylabel('Hits Success Rate (%)',fontsize=16)
plt.xlabel('Noise Scale (Arb.Units)', fontsize=16)
barplot.set_xticklabels((f'0.05\n n={len( df_filtered_112[df_filtered_112["Noise"]==0.05])}', 
                         f'0.1\n n={len(df_filtered_112[df_filtered_112["Noise"]==0.1])}', 
                         f'0.5\n n={len(df_filtered_112[df_filtered_112["Noise"]==0.5])}', 
                         f'1\n n={len(df_filtered_112[df_filtered_112["Noise"]==1])}'), fontsize=14)

# Show the plot
plt.ylim((0,0.07))
plt.tick_params(axis='y', labelsize=14)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/barplot_noisescale.png')


# Load data
df = pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/rmsd_noise_scale.csv')
print(df)
# Prepare data for plotting
df_filtered_112 = df[df['campaign'].str.contains('run_112_noise_scale')]
df_filtered_112['Noise'] = df_filtered_112['campaign'].str.extract('run_112_noise_scale_(\d+)').astype(str)
df_filtered_112.loc[df_filtered_112['Noise'] == '005', 'Noise'] = 0.05
df_filtered_112.loc[df_filtered_112['Noise'] == '01', 'Noise'] = 0.1
df_filtered_112.loc[df_filtered_112['Noise'] == '05', 'Noise'] = 0.5
df_filtered_112.loc[df_filtered_112['Noise'] == '1', 'Noise'] = 1.0












plt.figure(figsize=(12,10))
sns.set_theme(style='whitegrid')
sns.pointplot(data=df_filtered_112, x='Noise', y='RMSD', hue= 'Noise', palette='crest', s=100,capsize=0.2, ci=95 )
plt.title('RMSD vs Noise Scale', fontsize=20)
plt.xlabel('Noise Scale (Arb.U)', fontsize=16)
plt.ylabel(f'RMSD ($\AA$)', fontsize=16)
plt.legend(title='Noise Scale')
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/Noise_Scale_RMSD.png')

X=[0.05,0.1,0.5,1]
y=[]
for noise in df_filtered_112['Noise'].unique():
    y.append(np.mean(df_filtered_112['RMSD'][df_filtered_112['Noise']==noise]))

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)

print(' PEARSON CORRELATION COEFFICIENT VALUE:', r_value)

plt.figure(figsize=(8,10))
violinplot=sns.violinplot(data=df_filtered_112, x='Noise', y='RMSD', palette='crest')
annotator=Annotator(
    violinplot,
    data=df_filtered_112,
    x='Noise',
    y='RMSD',
    pairs=[
        (0.05,0.1),(0.05,0.5),(0.05,1),
        (0.1,0.5),(0.1,1),
        (0.5,1)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('RMSD vs Noise Scale', fontsize=20)
plt.ylabel(f'RMSD ($\AA$)', fontsize=16)
plt.xlabel('Noise Scale (Arb.U)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_rmsd_noisescale.png')

#df_filter hits
df_hits=df_filtered_112[(df_filtered_112['pae_interaction'] < 10) & (df_filtered_112['plddt_binder']>80)]
#A violinplot in case it is interesting 
plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_hits, x='Noise', y='pae_interaction', palette='crest')
annotator=Annotator(
    violinplot,
    data=df_hits,
    x='Noise',
    y='pae_interaction',
    pairs=[
        (0.05,0.10),(0.05,0.5),(0.05,1),
        (0.1,0.5),(0.1,1),
        (0.5,1)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Scale vs PAE interaction', fontsize=20)
plt.xlabel('Noise Scale (Arb.Units)', fontsize=16)
plt.ylabel(f'PAE interaction ($\AA$)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_pae_interaction_noise_scale.png')

plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_hits, x='Noise', y='plddt_binder', palette='crest')
annotator=Annotator(
    violinplot,
    data=df_hits,
    x='Noise',
    y='plddt_binder',
    pairs=[
        (0.05,0.10),(0.05,0.5),(0.05,1),
        (0.1,0.5),(0.1,1),
        (0.5,1)
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Scale vs pLDDT binder', fontsize=20)
plt.xlabel('Noise Scale (Arb.Units)',fontsize=16)
plt.ylabel(f'pLDDT binder (Arb.Units)', fontsize=16)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_plddt_binder_noise_scale.png')






df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv', header=0)

df_filtered_112 = df[df['campaign'].str.contains('run_112_noise_scale')]
df_filtered_112['Noise'] = df_filtered_112['campaign'].str.extract('run_112_noise_scale_(\d+)').astype(str)
df_filtered_112.loc[df_filtered_112['Noise'] == '005', 'Noise'] = 0.05
df_filtered_112.loc[df_filtered_112['Noise'] == '01', 'Noise'] = 0.1
df_filtered_112.loc[df_filtered_112['Noise'] == '05', 'Noise'] = 0.5
df_filtered_112.loc[df_filtered_112['Noise'] == '1', 'Noise'] = 1.0

# Original data point
original = {
    'pae_interaction': [7.383],
    'plddt_binder': [88.05]
}


df_original = pd.DataFrame(original)

# Set the style
sns.set(style="whitegrid")
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df_filtered_112,x='pae_interaction',y='plddt_binder', hue='Noise', s=100, palette='crest')
sns.scatterplot(data=df_original, x='pae_interaction', y='plddt_binder', color='red', s=100, edgecolors='k', label='Original')

# Add titles and labels
plt.title('PAE Interaction vs. pLDDT Binder Run_112', fontsize=20)
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
plt.legend(title='Noise Scale')

# Add grid
plt.grid(True)

# Show the plot
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseScaleRun112Figure.png')


#Figure for run_264

df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv', header=0)

df_filtered_264 = df[df['campaign'].str.contains('run_264_noise_scale')]
df_filtered_264['Noise'] = df_filtered_264['campaign'].str.extract('run_264_noise_scale_(\d+)').astype(str)
df_filtered_264.loc[df_filtered_264['Noise'] == '005', 'Noise'] = 0.05
df_filtered_264.loc[df_filtered_264['Noise'] == '01', 'Noise'] = 0.1
df_filtered_264.loc[df_filtered_264['Noise'] == '05', 'Noise'] = 0.5
df_filtered_264.loc[df_filtered_264['Noise'] == '1', 'Noise'] = 1.0

# Original data point
original = {
    'pae_interaction': [7.383],
    'plddt_binder': [88.05]
}

df_005ns=df_filtered_264[df_filtered_264['Noise']==0.05]
df_01ns=df_filtered_264[df_filtered_264['Noise']==0.1]
df_05ns=df_filtered_264[df_filtered_264['Noise']==0.5]
df_1ns=df_filtered_264[df_filtered_264['Noise']==1]

# Convert the original data to a DataFrame
df_original = pd.DataFrame(original)

# Set the style
sns.set(style="whitegrid")

# Create a scatter plot
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df_filtered_264, x='pae_interaction', y='plddt_binder', hue='Noise', s=100, palette='crest')
sns.scatterplot(data=df_original, x='pae_interaction', y='plddt_binder', color='red', s=100, edgecolors='k', label='Original')

# Add titles and labels
plt.title('PAE Interaction vs. pLDDT Binder Run_264', fontsize=20)
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
plt.legend()

# Add grid
plt.grid(True)

# Show the plot
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseScaleRun264Figure.png')

