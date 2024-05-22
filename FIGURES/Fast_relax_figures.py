import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

# Load the dataset
df = pd.read_csv('/emdata/cchacon/RFD_PD_test/output_af2.csv')

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
plt.text(3, 82, 'Hits Zone', color='black', fontsize=12)
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

percentage_FR=positives_FR/(positives_FR+negatives_FR)
percentage_noFR=positives_noFR/(positives_noFR+negatives_noFR)
# Prepare data for plotting
data = {'Campaign': ['Fast Relax', 'No Fast Relax'], 'Percentage': [percentage_FR, percentage_noFR]}

# Create a DataFrame
df_plot = pd.DataFrame(data)

# Set the style
sns.set(style="whitegrid")

# Create a bar plot
plt.figure(figsize=(8, 6))
barplot = sns.barplot(x='Campaign', y='Percentage', data=df_plot, palette=['blue', 'green'])

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
plt.gca().text(0.02, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
               verticalalignment='top', bbox=props)


# Add titles and labels
plt.title('Fast Relax vs No Fast Relax success rates')
plt.ylabel('Hits success rate')
# Show the plot
plt.show()