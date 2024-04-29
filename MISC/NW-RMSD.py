import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# Load data from CSV files
rmsd_df = pd.read_csv('rmsd.csv', header=0)
alignment_df = pd.read_csv('pairwise_alignment.csv', header=0)
metrics_df=pd.read_csv('../hits.csv', header=0)


# Merge DataFrames on 'query' column
merged_df_1 = pd.merge(rmsd_df, alignment_df, on='query')
merged_df=pd.merge(merged_df_1, metrics_df, on='query')
filtered_df = merged_df[merged_df['Global RMSD'] != 0]

# Create a scatter plot
plt.figure(figsize=(10, 6))  # Set figure size

# Scatter plot with labeled points
for query, group in filtered_df.groupby('query'):
    plt.scatter(group['Global RMSD'], group['score'], label=query)

# Add labels and title
plt.xlabel('Global RMSD')
plt.ylabel('Score')
plt.title('Score vs Global RMSD')

# Perform linear regression
X = filtered_df['Global RMSD'].values.reshape(-1, 1)
y = filtered_df['score'].values.reshape(-1, 1)
regression_model = LinearRegression().fit(X, y)

# Plot regression line
x_values = np.linspace(min(filtered_df['Global RMSD']), max(filtered_df['Global RMSD']), 100)
y_values = regression_model.predict(x_values.reshape(-1, 1))
plt.plot(x_values, y_values, color='red', linestyle='--', label='Regression Line')

# Calculate R-squared
r_squared = regression_model.score(X, y)
plt.text(0.80, 0.95, f'R-squared = {r_squared:.2f}', transform=plt.gca().transAxes, ha='left', va='top', fontsize=12)

# Add legend
plt.legend(title='Query', bbox_to_anchor=(1, 1), loc='upper left')

# Show grid
plt.grid(True)

# Show plot
plt.show()

#PLot for pae_interaction vs rmsd

plt.figure(figsize=(10, 6))  # Set figure size
for query, group in filtered_df.groupby('query'):
    plt.scatter(group['Global RMSD'], group['pae_interaction'], label=query)
plt.xlabel('Global RMSD')
plt.ylabel('par_interaction')
plt.title('Pae_interaction vs Global RMSD')
# Perform linear regression
X = filtered_df['Global RMSD'].values.reshape(-1, 1)
y = filtered_df['pae_interaction'].values.reshape(-1, 1)
regression_model = LinearRegression().fit(X, y)
slope = regression_model.coef_[0][0] 

print('The slope is: ', slope)
# Plot regression line
x_values = np.linspace(min(filtered_df['Global RMSD']), max(filtered_df['Global RMSD']), 100)
y_values = regression_model.predict(x_values.reshape(-1, 1))
plt.plot(x_values, y_values, color='red', linestyle='--', label='Regression Line')

# Calculate R-squared
r_squared = regression_model.score(X, y)
plt.text(0.80, 0.95, f'R-squared = {r_squared:.2f}', transform=plt.gca().transAxes, ha='left', va='top', fontsize=12)

# Add legend
plt.legend(title='Query', bbox_to_anchor=(1, 1), loc='upper left')

# Show grid
plt.grid(True)

# Show plot
plt.show()

#Plot for pae_interaction vs score

plt.figure(figsize=(10, 6))  # Set figure size
for query, group in filtered_df.groupby('query'):
    plt.scatter(group['score'], group['pae_interaction'], label=query)
plt.xlabel('score')
plt.ylabel('par_interaction')
plt.title('Pae_interaction vs score')
# Perform linear regression
X = filtered_df['score'].values.reshape(-1, 1)
y = filtered_df['pae_interaction'].values.reshape(-1, 1)
regression_model = LinearRegression().fit(X, y)
slope = regression_model.coef_[0][0] 

print('The slope is: ', slope)
# Plot regression line
x_values = np.linspace(min(filtered_df['score']), max(filtered_df['score']), 100)
y_values = regression_model.predict(x_values.reshape(-1, 1))
plt.plot(x_values, y_values, color='red', linestyle='--', label='Regression Line')

# Calculate R-squared
r_squared = regression_model.score(X, y)
plt.text(0.80, 0.95, f'R-squared = {r_squared:.2f}', transform=plt.gca().transAxes, ha='left', va='top', fontsize=12)

# Add legend
plt.legend(title='Query', bbox_to_anchor=(1, 1), loc='upper left')

# Show grid
plt.grid(True)

# Show plot
plt.show()


#plddt binder vs score 


plt.figure(figsize=(10, 6))  # Set figure size
for query, group in filtered_df.groupby('query'):
    plt.scatter(group['Global RMSD'], group['plddt_binder'], label=query)
plt.xlabel('Global RMSD')
plt.ylabel('plddt_binder')
plt.title('Plddt_binder vs Global RMSD')
# Perform linear regression
X = filtered_df['Global RMSD'].values.reshape(-1, 1)
y = filtered_df['plddt_binder'].values.reshape(-1, 1)
regression_model = LinearRegression().fit(X, y)
slope = regression_model.coef_[0][0] 

print('The slope is: ', slope)
# Plot regression line
x_values = np.linspace(min(filtered_df['Global RMSD']), max(filtered_df['Global RMSD']), 100)
y_values = regression_model.predict(x_values.reshape(-1, 1))
plt.plot(x_values, y_values, color='red', linestyle='--', label='Regression Line')

# Calculate R-squared
r_squared = regression_model.score(X, y)
plt.text(0.80, 0.95, f'R-squared = {r_squared:.2f}', transform=plt.gca().transAxes, ha='left', va='top', fontsize=12)

# Add legend
plt.legend(title='Query', bbox_to_anchor=(1, 1), loc='upper left')

# Show grid
plt.grid(True)

# Show plot
plt.show()