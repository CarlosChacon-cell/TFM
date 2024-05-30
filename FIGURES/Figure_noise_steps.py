import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import matplotlib.patches as patches
import scipy
from statannotations.Annotator import Annotator

# Load data
df = pd.read_csv('/emdata/cchacon/RFD_PD_test/pdbs/rmsd_noise_steps.csv')

# Prepare data for plotting
df_filtered = df[df['campaign'].str.contains('campaign_run_112_ns')]
df_filtered['Noise'] = df_filtered['campaign'].str.extract('campaign_run_112_ns(\d+)').astype(int)

# Initialize the plot
plt.figure(figsize=(10, 8))

# Plot with seaborn
sns.scatterplot(data=df_filtered, x='Noise', y='RMSD', hue='Noise', palette='viridis', s=100)

# Prepare data for regression with all the points
# X = df_filtered['Noise'].tolist()
# y = df_filtered['RMSD'].tolist()
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
# r_squared = r_value ** 2
# y_pred = [intercept + slope * value for value in X]

X=[5,10,20,30]
y=[]
for noise in df_filtered['Noise'].unique():
    y.append(np.mean(df_filtered['RMSD'][df_filtered['Noise']==noise]))

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
y_pred = [intercept + slope * value for value in X]

# Plot the regression line
sns.lineplot(x=X, y=y_pred, label=f'$r$={r_value:.2f}', color='red')

# Add title and labels
plt.title('RMSD vs. Noise Levels with Linear Fits and $R^2$ Values', fontsize=16)
plt.xlabel('Noise Level', fontsize=14)
plt.ylabel('RMSD', fontsize=14)
plt.legend(title='Regression Line')
plt.grid(True)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseStepsRMSD.png')

plt.figure(figsize=(8,10))
boxplot=sns.boxplot(data=df_filtered, x='Noise', y='RMSD', palette='viridis')
annotator=Annotator(
    boxplot,
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
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/boxplot_rmsd_noisesteps.png')

#df_filter hits
df_hits=df_filtered[(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder']>80)]
#A boxplot in case it is interesting 
plt.figure(figsize=(10,12))
boxplot=sns.boxplot(data=df_hits, x='Noise', y='pae_interaction', palette='viridis')
annotator=Annotator(
    boxplot,
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

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/boxplot_pae_interaction_noise_steps.png')




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
plt.scatter(df_5ns['pae_interaction'], df_5ns['plddt_binder'], color='blue', label='5 Noise Steps')
plt.scatter(df_10ns['pae_interaction'], df_10ns['plddt_binder'], color='purple', label='10 Noise Steps')
plt.scatter(df_20ns['pae_interaction'], df_20ns['plddt_binder'], color='orange', label='20 Noise Steps')
plt.scatter(df_30ns['pae_interaction'], df_30ns['plddt_binder'], color='green', label='30 Noise Steps')
plt.scatter(df_original['pae_interaction'], df_original['plddt_binder'], color='red', s=100, label='Original', edgecolors='black')

# Add titles and labels
plt.title('PAE Interaction vs. pLDDT Binder')
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
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseStepsRun264Figure.png')

