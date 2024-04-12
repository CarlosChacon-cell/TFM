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


for campaign in set(df['campaign']):
    total_designs=len(df[df['campaign']==campaign])
    success_designs=len(filtered_df[filtered_df['campaign']==campaign])
    success_rate=success_designs/total_designs

    print(campaign,' Total designs: ',total_designs, ' Success: ', success_designs, ' Success_rate: ',success_rate*100 )