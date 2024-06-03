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

X=[5,10,20,30]
y=[]
for noise in df_filtered['Noise'].unique():
    y.append(np.mean(df_filtered['RMSD'][df_filtered['Noise']==noise]))

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
y_pred = [intercept + slope * value for value in X]

# Initialize the plot

plt.figure(figsize=(12,10))
sns.set_theme(style='whitegrid')
sns.pointplot(data=df_filtered, x='Noise', y='RMSD', hue= 'Noise', palette='crest', s=100,capsize=0.2, ci=95)


# Plot the regression line
# sns.lineplot(x=X, y=y_pred, label=f'$r$={r_value:.2f}', color='red', ax=ax) ##Doesn't work##
plt.title('RMSD vs Noise Steps')
plt.xlabel('Noise Steps (#)')
plt.ylabel(f'RMSD ($\AA$)')
plt.legend(title='Noise Steps')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseStepsRMSD.png')

plt.figure(figsize=(8,10))
violinplot=sns.violinplot(data=df_filtered, x='Noise', y='RMSD', palette='viridis')
annotator=Annotator(
    violinplot,
    data=df_filtered,
    x='Noise',
    y='RMSD',
    pairs=[
        (5,10),(5,20),(5,30),
        (10,20),(10,30),
        (20,30)
    ]
)
annotator.configure(
    test='t-test_ind',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('RMSD vs Noise Steps')
plt.ylabel(f'RMSD ($\AA$)')
plt.xlabel('Noise Steps(#)')
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_rmsd_noisesteps.png')

#df_filter hits
df_hits=df_filtered[(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder']>80)]
d_stat,p_value=ks_2samp(df_hits['RMSD'][df_hits['Noise']==10], df_hits['RMSD'][df_hits['Noise']==5])
print('PVALUE BY KS IS: ',  p_value)

#A violinplot in case it is interesting 
plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_hits, x='Noise', y='pae_interaction', palette='viridis')
annotator=Annotator(
    violinplot,
    data=df_hits,
    x='Noise',
    y='pae_interaction',
    pairs=[
        (5,10),(5,20),(5,30),
        (10,20),(10,30),
        (20,30)
    ]
)
annotator.configure(
    test='t-test_ind',
    text_format='star',
    loc='inside',
    comparisons_correction="bonferroni")
    
annotator.apply_and_annotate()
plt.title('Noise Steps vs Pae_interaction')
plt.xlabel('Noise Steps (#)')
plt.ylabel(f'Pae_interaction ($\AA$)')

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_pae_interaction_noise_steps.png')

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
plt.title('PAE Interaction vs. pLDDT Binder run_264')
plt.xlabel('PAE Interaction')
plt.ylabel('pLDDT Binder')

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

