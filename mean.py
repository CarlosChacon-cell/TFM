import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('merged_output_pd.csv')

# Group the DataFrame by the 'descriptor' column and calculate the mean for 'plddt_binder' and 'pae_interaction'
mean_values = df.groupby('campaign')[['plddt_binder', 'pae_interaction']].std()

print(mean_values)
