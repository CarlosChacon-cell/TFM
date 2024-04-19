import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('/emdata/cchacon/TRF1/20240415_campaign_hotspotDIM_TRF1/pae_local_global.csv')

global_pae = df['pae_interaction_global']
local_pae = df['pae_interaction_local']

# Extract unique protein names and assign a color to each protein
unique_proteins = df['protein_name'].unique()
colors = plt.cm.get_cmap('tab20', len(unique_proteins))  # Using a colormap for distinct colors

plt.figure(figsize=(8, 6))

# Scatter plot for global PAE with colors based on protein names
plt.scatter(np.zeros(len(global_pae)), global_pae, label='Global Pae', color=colors(df['protein_name'].apply(lambda x: np.where(unique_proteins == x)[0][0])), s=100)  # Increased point size (s=100)

# Scatter plot for local PAE with colors based on protein names
plt.scatter(np.ones(len(local_pae)), local_pae, label='Local Pae', color=colors(df['protein_name'].apply(lambda x: np.where(unique_proteins == x)[0][0])), s=100)  # Increased point size (s=100)

for i in range(len(global_pae)):
    plt.plot([0, 1], [global_pae[i], local_pae[i]], c='k', lw=1)  # Thinner line (lw=1)

plt.xticks([0, 1], ['Global PAE', 'Local PAE'])
plt.ylabel('pae_interaction')
plt.title('Global vs local PAE interaction')
plt.grid(True)
plt.show()



df2=pd.read_csv('/emdata/cchacon/TRF1/20240415_campaign_hotspotDIM_TRF1/close_residues.csv')

length=df2['length']
interacting_residues=df2['interacting_residues']

interacting_surface=[]
for i in range(len(length)):
    interacting_surface.append(len(interacting_residues)/int(length))

    