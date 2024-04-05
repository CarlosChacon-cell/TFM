import pandas as pd


# Load the CSV file into a DataFrame
df = pd.read_csv('output_af2.csv')

# Filter the DataFrame based on conditions
filtered_df = df[(df['pae_interaction'] < 10) & (df['plddt_binder'] > 80)]

# Group the filtered DataFrame by the 'descriptor' column and calculate the mean for 'plddt_binder' and 'pae_interaction'
mean_values = filtered_df.groupby('campaign')[['plddt_binder', 'pae_interaction']].mean()

print('mean values: \n', mean_values)

median_values = filtered_df.groupby('campaign')[['plddt_binder', 'pae_interaction']].median()

print('median values: \n', median_values)

total_designs=df.groupby('campaign').len()
success_designs=filtered_df.groupby('campaign').len()

success_rate=success_designs/total_designs

print('success rates: \n', success_rate)