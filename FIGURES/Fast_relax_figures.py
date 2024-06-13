import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from statannotations.Annotator import Annotator
import numpy as np
from scipy.optimize import curve_fit



# Load the dataset
df = pd.read_csv('/emdata/cchacon/RFD_PD_test/output_af2.csv')
df_filtered=df[df['campaign'].str.contains('FR')]
df_filtered=df_filtered[(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder'] > 80)]


# Filter the data
df_FR = df[df['campaign'] == 'campaign_run_112_FR']
df_noFR = df[df['campaign'] == 'campaign_run_112_noFR']

# Original data point
original = {
    'pae_interaction': [8.97],
    'plddt_binder': [81.822]
}

# Convert the original data to a DataFrame
df_original = pd.DataFrame(original)

# Set the style
sns.set(style="whitegrid")

# Create a scatter plot
plt.figure(figsize=(10, 8))
plt.scatter(df_FR['pae_interaction'], df_FR['plddt_binder'], color='blue', label='Fast Relax')
plt.scatter(df_noFR['pae_interaction'], df_noFR['plddt_binder'], color='green', label='No Fast Relax')
plt.scatter(df_original['pae_interaction'], df_original['plddt_binder'], color='red', s=100, label='Original', edgecolors='black')

# Add titles and labels
plt.title('PAE Interaction vs. pLDDT Binder')
plt.xlabel('PAE Interaction')
plt.ylabel('pLDDT Binder')

# Highlight specific region with a square
rect = patches.Rectangle((0, 80), 10, 10, linewidth=2, edgecolor='black', facecolor='none')
plt.gca().add_patch(rect)
plt.text(4, 81, 'Hits Zone', color='black', fontsize=12)
# Invert the x-axis
plt.gca().invert_xaxis()

# Add legend
plt.legend()

# Add grid
plt.grid(True)

# Show the plot
# plt.show()

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/FastRelaxFigure.png')

positives_FR=len(df_FR[(df_FR['plddt_binder'] > 80) & (df_FR['pae_interaction'] < 10)])
negatives_FR = len(df_FR) - positives_FR

positives_noFR=len(df_noFR[(df_noFR['plddt_binder'] > 80) & (df_noFR['pae_interaction'] < 10)])
negatives_noFR = len(df_noFR) - positives_noFR

percentage_FR=positives_FR/(positives_FR+negatives_FR) *100
percentage_noFR=positives_noFR/(positives_noFR+negatives_noFR) *100
# Prepare data for plotting
data = {'Campaign': ['Fast Relax', 'No Fast Relax'], 'Percentage': [percentage_FR, percentage_noFR]}

# Create a DataFramea
df_plot = pd.DataFrame(data)

# Set the style
sns.set(style="whitegrid")

# Create a bar plot
plt.figure(figsize=(10,12))
barplot = sns.barplot(x='Campaign', y='Percentage', data=df_plot, palette='crest')

# Add significance bracket with an asterisk
x1, x2 = 0, 1   # x-coordinates of the bars
y, h, col = max(df_plot['Percentage']) + 0.01, 0.02, 'k'  # y-coordinate, height of the bracket, color
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
# Add a box with the p-value in the upper left corner
p_value = 0.048  # Replace with the actual p-value when available
textstr = f'p-value = {p_value}'
props = dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white')

# Place the box in the upper left corner
# plt.gca().text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
#                verticalalignment='top', bbox=props)


# Add titles and labels
plt.title('Fast Relax vs No Fast Relax Success Rates', fontsize=26)
plt.ylabel('Hits Success Rate (%)', fontsize=22)
plt.xlabel('Fast Relax Protocol', fontsize=22)
barplot.set_xticklabels((f'Fast Relax\n n={positives_FR + negatives_FR}', f'No Fast Relax\n n={positives_noFR+negatives_noFR}'), fontsize=20)
# Show the plot
plt.tick_params(axis='y', labelsize=20)

plt.ylim((0,7))
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/FastRelaxvsNoFRSuccessRates.png')

plt.figure(figsize=(10,12))
hits_FR=len(df_filtered[df_filtered['campaign']=='campaign_run_112_FR'])
hits_noFR=len(df_filtered[df_filtered['campaign']=='campaign_run_112_noFR'])
violinplot=sns.violinplot(data=df_filtered, x='campaign', y='pae_interaction', palette='crest')
violinplot.set_xticklabels((f'Fast Relax\n n={hits_FR}', f'No Fast Relax\n n={hits_noFR}'), fontsize=20)
annotator=Annotator(
    violinplot,
    data=df_filtered,
    x='campaign',
    y='pae_interaction',
    pairs=[
        ('campaign_run_112_FR', 'campaign_run_112_noFR')
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="holm-bonferroni",
    fontsize=18)
    
annotator.apply_and_annotate()
plt.title('PAE vs Fast Relax', fontsize=26)
plt.ylabel(f'PAE interaction ($\AA$)', fontsize=22)
plt.xlabel('Fast Relax Protocol', fontsize=22)

plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)

plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_pae_interaction_FR.png')

plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_filtered, x='campaign', y='plddt_binder', palette='crest')
violinplot.set_xticklabels((f'Fast Relax\n n={hits_FR}', f'No Fast Relax\n n={hits_noFR}'), fontsize=20)
annotator=Annotator(
    violinplot,
    data=df_filtered,
    x='campaign',
    y='plddt_binder',
    pairs=[
        ('campaign_run_112_FR', 'campaign_run_112_noFR')
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="holm-bonferroni",
    fontsize=18)
    
annotator.apply_and_annotate()

plt.title('pLDDT vs _Fast Relax', fontsize=26)
plt.ylabel(f'pLDDT binder', fontsize=22)
plt.xlabel('Fast Relax Protocol', fontsize=22)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_plddt_binder_FR.png')

rmsd_fr_path='/emdata/cchacon/RFD_PD_test/campaign_run_112_FR/hits/rmsd.csv'
rmsd_nofr_path='/emdata/cchacon/RFD_PD_test/campaign_run_112_noFR/hits/rmsd.csv'

df_rmsd_fr=pd.read_csv(rmsd_fr_path, header=0)
df_rmsd_fr['FR']='Fast Relax'
df_rmsd_nofr=pd.read_csv(rmsd_nofr_path, header=0)
df_rmsd_nofr['FR']='No Fast Relax'

df_rmsd=pd.concat([df_rmsd_fr, df_rmsd_nofr])
plt.figure(figsize=(10,12))
violinplot=sns.violinplot(data=df_rmsd, x='FR', y='Global RMSD', palette='crest')
violinplot.set_xticklabels((f'Fast Relax\n n={hits_FR}', f'No Fast Relax\n n={hits_noFR}'), fontsize=20)

annotator=Annotator(
    violinplot,
    data=df_rmsd,
    x='FR',
    y='Global RMSD',
    pairs=[
        ('Fast Relax', 'No Fast Relax')
    ]
)
annotator.configure(
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    comparisons_correction="holm-bonferroni",
    fontsize=18)
    
annotator.apply_and_annotate()
plt.title('RMSD vs _Fast Relax', fontsize=26)
plt.ylabel(f'RMSD ($\AA$)', fontsize=22)
plt.xlabel('Fast Relax Protocol', fontsize=22)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/violinplot_rmsd_FR.png')