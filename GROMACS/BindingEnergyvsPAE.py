import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import re

# Load the CSV file
df = pd.read_csv('deltavspae.csv', header=None, index_col=None)
df = df.rename(columns={0: 'protein_name', 1: 'Delta', 2: 'Pae_interaction'})
df_local=pd.read_csv('pae_local_global.csv', header=0)
merged_df=pd.merge(df,df_local, on='protein_name' , how='left')
print(merged_df)

# Create a colormap with as many colors as there are unique proteins
colors = plt.cm.get_cmap('tab20', len(df['protein_name'].unique()))
pattern = r'run_\d+_design_\d+(?:_PD)?'
# Scatter plot
plt.figure(figsize=(8, 6))
for i, protein in enumerate(merged_df['protein_name'].unique()):
    subset_df = merged_df[merged_df['protein_name'] == protein]
    protein_name=re.search(pattern, protein).group()
    plt.scatter(subset_df['Delta'], subset_df['Pae_interaction'], label=protein_name, color=colors(i))

# Fit a line (linear regression)
X = merged_df[['Delta']]
y = merged_df['Pae_interaction']
model = LinearRegression().fit(X, y)
y_pred = model.predict(X)
r_squared = r2_score(y, y_pred)
# Plot the regression line
plt.plot(X, y_pred, color='r', linewidth=1, linestyle='--', label=f'R² = {round(r_squared,2)}')

# Calculate R-squared


plt.xlabel('Delta')
plt.ylabel('Pae_interaction')
plt.title('Binding Energy vs Pae_interaction')
plt.grid(True)
plt.legend(loc='upper left')
plt.savefig('DeltaVsPae.png')

plt.figure(figsize=(8, 6))
for i, protein in enumerate(merged_df['protein_name'].unique()):
    subset_df = merged_df[merged_df['protein_name'] == protein]
    protein_name=re.search(pattern, protein).group()
    plt.scatter(subset_df['Delta'], subset_df['pae_interaction_local'], label=protein_name, color=colors(i))

# Fit a line (linear regression)
X = merged_df[['Delta']]
y = merged_df['pae_interaction_local']
model = LinearRegression().fit(X, y)
y_pred = model.predict(X)
r_squared = r2_score(y, y_pred)
# Plot the regression line
plt.plot(X, y_pred, color='r', linewidth=1, linestyle='--', label=f'R² = {round(r_squared,2)}')

# Calculate R-squared


plt.xlabel('Delta')
plt.ylabel('Pae_interaction_local')
plt.title('Binding Energy vs Pae_interaction_local')
plt.grid(True)
plt.legend(loc='upper left')
plt.savefig('DeltaVsPae_local.png')
