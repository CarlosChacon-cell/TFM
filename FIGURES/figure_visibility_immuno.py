import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
df_1 = pd.read_csv('/emdata/cchacon/inmuno/Scaffold_run_91/run_5/Visibility.csv')
df_2 = pd.read_csv('/emdata/cchacon/inmuno/Scaffold_run_91/run_2/Visibility.csv')
df_3 = pd.read_csv('/emdata/cchacon/inmuno/Scaffold_run_91/Visibility.csv')

# Add identifier for original data
df_3['file'] = 'Original'

# Concatenate dataframes
df_concat = pd.concat([df_1, df_2, df_3])

# Create the plot
plt.figure(figsize=(12, 8))
sns.histplot(data=df_concat, x='Visibility', palette='Set1', element='step', stat='density')

# Add vertical lines for the 'Original' data
original_visibility = df_concat['Visibility'][df_concat['file'] == 'Original']
plt.vlines(x=original_visibility, ymin=0, ymax=plt.gca().get_ylim()[1], color='red', linestyles='--', label='Original Visibility')

# Add a title and labels
plt.title('Visibility Distribution', fontsize=26)
plt.xlabel('Visibility', fontsize=22)
plt.ylabel('Density', fontsize=22)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(axis='x', labelsize=20)


# Adjust the legend
plt.legend(fontsize=14)

# Show the plot
plt.savefig('/home/cchacon/Carlos_scripts/FIGURES/VisibilityDistributionScaffold94.png')
