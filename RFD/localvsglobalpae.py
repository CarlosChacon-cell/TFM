import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns 

file_path=glob.glob('pae_local_global.csv')[0]
print(file_path)

df = pd.read_csv(file_path)

global_pae = df['pae_interaction_global']
local_pae = df['pae_interaction_local']

# Extract unique protein names and assign a color to each protein
unique_proteins = df['protein_name'].unique()
colors = plt.cm.get_cmap('tab20', len(unique_proteins))  # Using a colormap for distinct colors

# plt.figure(figsize=(8, 6))

# # Scatter plot for global PAE with colors based on protein names
# plt.scatter(np.zeros(len(global_pae)), global_pae, label='Global Pae', color=colors(df['protein_name'].apply(lambda x: np.where(unique_proteins == x)[0][0])), s=100)  # Increased point size (s=100)

# # Scatter plot for local PAE with colors based on protein names
# plt.scatter(np.ones(len(local_pae)), local_pae, label='Local Pae', color=colors(df['protein_name'].apply(lambda x: np.where(unique_proteins == x)[0][0])), s=100)  # Increased point size (s=100)

# for i in range(len(global_pae)):
#     plt.plot([0, 1], [global_pae[i], local_pae[i]], c='k', lw=1)  # Thinner line (lw=1)

# plt.xticks([0, 1], ['Global PAE', 'Local PAE'])
# plt.ylabel('pae_interaction')
# plt.title('Global vs local PAE interaction')
# plt.grid(True)
# plt.show()



# df2=pd.read_csv('/emdata/cchacon/TRF1/20240415_campaign_hotspotDIM_TRF1/close_residues.csv')

# length=df2['length'] 
# for i in range(len(length)):
#     interacting_surface.append(len(interacting_residues)/int(length))
counter=0
total=0
differences=[]
for i in range(len(df['protein_name'])):
    if df['pae_interaction_local'][i] > df['pae_interaction_global'][i]:
        differences.append(df['pae_interaction_local'][i]-df['pae_interaction_global'][i])
        counter+=1
        total+=1
    else:
        differences.append(df['pae_interaction_local'][i]-df['pae_interaction_global'][i])
        total+=1

improve=100-(counter/total)*100 

print('Local pae_interaction predicts an improvement for: ', improve, ' of the designs')
print('')
print('Local pae_interaction - global_pae_interaction: ', np.mean(differences))



print('Mean interacting surface: ', np.mean(df['interacting_surface']))


# Initialize the plot
plt.figure(figsize=(10, 8))

# Plot with seaborn
sns.scatterplot(data=df, y='pae_interaction_local', x='pae_interaction_global', hue='interacting_surface', palette='viridis', s=100)
plt.xlabel('Pae_interaction')
plt.ylabel('CUTRE')
plt.title('Pae_interaction vs CUTRE')
plt.grid(True)
plt.xlim((0,35))
plt.ylim((0,35))

plt.show()