import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import re

files_list=glob.glob('*/secondary_struct*')
pattern=r'.*/secondary_structure_(mut\d+)_(\d+).txt'
# Load the data
for file in files_list:
    file_path = file
    print(file_path)
    protein_name=re.search(pattern,file_path).group(1)
    
    dataset= re.search(pattern,file_path).group(2)

    df_data = pd.read_csv(file_path, header=None, sep='\s+', skiprows=range(10), names=['Residue_number', 'Uk', 'Emptyval', 'Frame', 'Type'])

    # Map secondary structure types to integers
    structure_mapping = {structure: idx for idx, structure in enumerate(sorted(df_data['Type'].unique()))}
    print(structure_mapping)
    df_data['Type_mapped'] = df_data['Type'].map(structure_mapping)

    # Pivot the dataframe to get the desired format for heatmap
    df_pivot_mapped = df_data.pivot(index='Residue_number', columns='Frame', values='Type_mapped')
    print(df_pivot_mapped)

    # Create a custom color palette
    palette = sns.color_palette("Set2", len(structure_mapping))

    # Create the heatmap without lines
    plt.figure(figsize=(15, 10))
    heatmap = sns.heatmap(df_pivot_mapped, cmap=palette, cbar_kws={'label': 'Secondary Structure Type'}, linewidths=0)

    # Create custom colorbar labels
    colorbar = heatmap.collections[0].colorbar
    colorbar.set_ticks(list(structure_mapping.values()))
    colorbar.set_ticklabels(list(structure_mapping.keys()))

    # Set labels and title
    plt.xlabel('Time (ps)')
    plt.ylabel('Residue Number')
    plt.title(f'Secondary Structure Heatmap of {protein_name} and dataset {dataset} ')

    # Show the plot
    plt.show()
