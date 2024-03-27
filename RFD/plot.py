

import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
filename = "merged_output_pd.csv"
data = pd.read_csv(filename)

# Extract data for plotting
pae_interaction = data['pae_interaction']
plddt_binder = data['plddt_binder']
line_filename = data['description']

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(pae_interaction, plddt_binder)

# Label each point with line filename
for i, txt in enumerate(line_filename):
    plt.annotate(txt.split('_design')[0], (pae_interaction[i], plddt_binder[i]))

# Set labels and title
plt.xlabel('PAE Interaction')
plt.ylabel('PLDDT Binder')
plt.title('PAE Interaction vs PLDDT Binder')

# Show plot
plt.grid(True)
plt.show()
