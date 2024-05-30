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

# Initialize the plot
plt.figure(figsize=(10, 8))

# Plot with seaborn
sns.scatterplot(data=df_filtered_112, x='Noise', y='RMSD', hue='Noise', palette='viridis', s=100)

# Prepare data for regression with all the points
# X = df_filtered_112['Noise'].tolist()
# y = df_filtered_112['RMSD'].tolist()
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
# r_squared = r_value ** 2
# y_pred = [intercept + slope * value for value in X]

X=[0.05,0.1, 0.5,1]
y=[]
for noise in df_filtered_112['Noise'].unique():
    y.append(np.mean(df_filtered_112['RMSD'][df_filtered_112['Noise']==noise]))

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, y)
y_pred = [intercept + slope * value for value in X]

# Plot the regression line
sns.lineplot(x=X, y=y_pred, label=f'$r$={r_value:.2f}', color='red')
print(r_value)
# Add title and labels
plt.title('RMSD vs. Noise Scale Levels', fontsize=16)
plt.xlabel('Noise Level', fontsize=14)
plt.ylabel('RMSD', fontsize=14)
plt.legend(title='Regression Line', loc='lower right')
plt.grid(True)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/Noise_Scale_RMSD.png')

plt.figure(figsize=(8,10))
boxplot=sns.boxplot(data=df_filtered_112, x='Noise', y='RMSD', palette='viridis')
annotator=Annotator(
    boxplot,
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
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/boxplot_rmsd_noisescale.png')




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
plt.scatter(df_005ns['pae_interaction'], df_005ns['plddt_binder'], color='blue', label='0.05 Noise Scale')
plt.scatter(df_01ns['pae_interaction'], df_01ns['plddt_binder'], color='purple', label='0.1 Noise Scale')
plt.scatter(df_05ns['pae_interaction'], df_05ns['plddt_binder'], color='orange', label='0.5 Noise Scale')
plt.scatter(df_1ns['pae_interaction'], df_1ns['plddt_binder'], color='green', label='1 Noise Scale')
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
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/NoiseScaleRun264Figure.png')

