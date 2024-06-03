import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
file_path = 'path_to_your_file/secondary_structure.txt'
df_data = pd.read_csv(file_path, header=None, sep='\s+', skiprows=range(10), names=['Residue_number', 'Uk', 'Emptyval', 'Frame', 'Type'])

# Map secondary structure types to integers
structure_mapping = {structure: idx for idx, structure in enumerate(df_data['Type'].unique())}
df_data['Type_mapped'] = df_data['Type'].map(structure_mapping)

# Pivot the dataframe to get the desired format for heatmap
df_pivot_mapped = df_data.pivot(index='Frame', columns='Residue_number', values='Type_mapped')

# Create a custom color palette
palette = sns.color_palette("tab10", len(structure_mapping))

# Create the heatmap without lines
plt.figure(figsize=(15, 10))
heatmap = sns.heatmap(df_pivot_mapped, cmap=palette, cbar_kws={'label': 'Secondary Structure Type'}, linewidths=0)

# Create custom colorbar labels
colorbar = heatmap.collections[0].colorbar
colorbar.set_ticks(list(structure_mapping.values()))
colorbar.set_ticklabels(list(structure_mapping.keys()))

# Set labels and title
plt.xlabel('Residue Number')
plt.ylabel('Frame')
plt.title('Secondary Structure Heatmap')

# Show the plot
plt.show()
