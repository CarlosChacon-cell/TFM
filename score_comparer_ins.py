import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

ins_scores = pd.read_csv('/emdata/cchacon/MD_RFD/instability_index/instability_index.csv', header=0, sep=',')
print(ins_scores)

af2_scores = pd.read_csv('/emdata/cchacon/NR_RFL_RFD/af2_scores.csv', header=0, sep=',')
print(af2_scores)

pattern = r'^run_\d+'

for i in range(len(af2_scores['description'])):
    descriptor = af2_scores['description'][i]
    target_name = re.search(pattern, descriptor)
    target_name = target_name.group(0)
    af2_scores['description'][i] = target_name

df = pd.merge(ins_scores, af2_scores, on='description')
print(df)

ax = df.plot.scatter(x='ins_index', y='plddt_binder')

# Perform linear regression
fit = np.polyfit(df['ins_index'], df['plddt_binder'], 1)
fit_fn = np.poly1d(fit)

# Plot the linear regression line
ax.plot(df['ins_index'], fit_fn(df['ins_index']), color='red', linestyle='dashed')

# Annotate each point with the target_name
for idx, row in df.iterrows():
    ax.annotate(row['description'], (row['ins_index'], row['plddt_binder']))

# Add labels and title
ax.set_xlabel('Ins Index')
ax.set_ylabel('PLDDT Binder')
ax.set_title('Scatter Plot with Linear Regression')

plt.show()
