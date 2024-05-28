import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns 

# Load the CSV file
file_path = glob.glob('pae_local_global.csv')[0]
print(file_path)

df = pd.read_csv(file_path)

# Extract global and local PAE interactions
global_pae = df['pae_interaction_global']
local_pae = df['pae_interaction_local']

# Calculate the differences and improvements
counter = 0
total = 0
differences = []

for i in range(len(df['protein_name'])):
    if df['pae_interaction_local'][i] >= df['pae_interaction_global'][i]:
        differences.append(df['pae_interaction_local'][i] - df['pae_interaction_global'][i])
        counter += 1
        total += 1
    else:
        differences.append(df['pae_interaction_local'][i] - df['pae_interaction_global'][i])
        total += 1

improve = 100 - (counter / total) * 100 

print('Local pae_interaction predicts an improvement for: ', improve, ' of the designs')
print('')
print('Local pae_interaction - global_pae_interaction: ', np.mean(differences))
print('Mean interacting surface: ', np.mean(df['interacting_surface']))

# Initialize the plot
plt.figure(figsize=(10, 8))

# Plot with seaborn
sns.scatterplot(data=df, y='pae_interaction_local', x='pae_interaction_global', hue='interacting_surface', palette='viridis', s=100)
plt.xlabel('Global PAE Interaction')
plt.ylabel('Local PAE Interaction')
plt.title('Global PAE Interaction vs Local PAE Interaction')
plt.grid(True)
plt.xlim((0, 30))
plt.ylim((0, 30))

# Add vertical and horizontal lines
plt.axvline(x=10, color='red', linestyle='--', label='Global threshold')
plt.axhline(y=10, color='blue', linestyle='--', label='CUTRE threshold')
plt.legend()
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/Figure_CUTREvsPaeInteraction.png')

plt.figure(figsize=(10,8))
df['Improvement']=df['pae_interaction_global']-df['pae_interaction_local']
sns.scatterplot(data=df, x='interacting_surface', y='Improvement', hue='pae_interaction_global', s=100)
plt.show()