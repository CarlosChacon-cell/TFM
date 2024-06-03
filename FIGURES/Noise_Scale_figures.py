import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import scipy
import matplotlib.patches as patches
from statannotations.Annotator import Annotator
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
plt.title('RMSD vs Noise Scale')
plt.xlabel('Noise Scale (Arb.U)')
plt.ylabel(f'RMSD ($\AA$)')
plt.legend(title='Noise Scale')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/Noise_Scale_RMSD.png')


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
    test='t-test_ind',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('RMSD vs Noise Scale')
plt.ylabel(f'RMSD ($\AA$)')
plt.xlabel('Noise Scale (Arb.U)')
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
    test='t-test_ind',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Scale vs Pae_interaction')
plt.xlabel('Noise Scale (Arb.Units)')
plt.ylabel(f'Pae_interaction ($\AA$)')

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_pae_interaction_noise_scale.png')

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
plt.title('PAE Interaction vs. pLDDT Binder Run_112')
plt.xlabel('PAE Interaction')
plt.ylabel('pLDDT Binder')

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
plt.title('PAE Interaction vs. pLDDT Binder Run_264')
plt.xlabel('PAE Interaction')
plt.ylabel('pLDDT Binder')

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

