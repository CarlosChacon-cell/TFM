import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Load the data
df = pd.read_csv('/emdata/cchacon/TRF1/20240412_campaign_hotspotA_B_myb/all_models.csv')

# Prepare the data
x = df['length'].values.reshape(-1, 1)
y = df['plddt_binder'].values

# Fit linear regression model
model = LinearRegression()
model.fit(x, y)

# Get the slope (m) and intercept (b) of the line
m = model.coef_[0]
b = model.intercept_

# Make predictions
y_pred = model.predict(x)

# Calculate R-squared
r_squared = r2_score(y, y_pred)

# Plotting
plt.figure(figsize=(8, 6))
plt.scatter(x, y, label='Original data')
plt.plot(x, y_pred, color='red', label=f'Linear Regression\ny = {m:.2f}x + {b:.2f}\nR-squared = {r_squared:.2f}')
plt.xlabel('length')
plt.ylabel('plddt_binder')
plt.title('Linear Regression')
plt.legend()
plt.grid(True)
plt.show()
