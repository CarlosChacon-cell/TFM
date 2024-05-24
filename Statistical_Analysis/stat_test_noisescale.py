from scipy.stats import chi2_contingency
import pandas as pd 
import numpy as np
from scipy.stats import ttest_ind
import argparse


df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv')
# Prepare data for plotting
df_filtered = df[df['campaign'].str.contains('run_112_noise_scale')]
df_filtered['Noise'] = df_filtered['campaign'].str.extract('run_112_noise_scale_(\d+)').astype(str)
print(df)

pae_interaction_ns005=df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['Noise']=='005')]
pae_interaction_ns01=df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['Noise']=='01')]
pae_interaction_ns05=df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['Noise']=='05')]
pae_interaction_ns1=df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['Noise']=='1')]

bonferroni_correction=0.05/6

# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns01, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 01 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 05 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns1, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns01, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 01 and 05 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns01, pae_interaction_ns1, equal_var=False)
print(f'T-Test for pae_interaction between 01 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns1, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 05 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')
if p_value < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')

positives_005ns=len(df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder'] > 80) & (df_filtered['Noise']=='005')])
total_005ns=len(df_filtered['pae_interaction'][(df_filtered['Noise']=='005')])

positives_01ns=len(df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder'] > 80) & (df_filtered['Noise']=='01')])
total_01ns=len(df_filtered['pae_interaction'][(df_filtered['Noise']=='01')])

positives_05ns=len(df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder'] > 80) & (df_filtered['Noise']=='05')])
total_05ns=len(df_filtered['pae_interaction'][(df_filtered['Noise']=='05')])

positives_1ns=len(df_filtered['pae_interaction'][(df_filtered['pae_interaction'] < 10) & (df_filtered['plddt_binder'] > 80) & (df_filtered['Noise']=='1')])
total_1ns=len(df_filtered['pae_interaction'][(df_filtered['Noise']=='1')])


# Create a contingency table
observed = [[positives_005ns, total_005ns - positives_005ns], [positives_01ns, total_01ns- positives_01ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 005 and 01 ns:", chi2)
print("P-value between 005 and 01 ns:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_05ns, total_05ns - positives_05ns], [positives_005ns, total_005ns- positives_005ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 005 and 05:", chi2)
print("P-value between 005 and 05:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_1ns, total_1ns - positives_1ns], [positives_005ns, total_005ns- positives_005ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 005 and 1:", chi2)
print("P-value between 005 and 1:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_1ns, total_1ns - positives_1ns], [positives_01ns, total_01ns- positives_01ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 01 and 1:", chi2)
print("P-value between 01 and 1:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')

# Create a contingency table
observed = [[positives_1ns, total_1ns - positives_1ns], [positives_05ns, total_05ns- positives_05ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 05 and 1:", chi2)
print("P-value between 05 and 1:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')

# Create a contingency table
observed = [[positives_05ns, total_05ns - positives_05ns], [positives_01ns, total_01ns- positives_01ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 01 and 05:", chi2)
print("P-value between 01 and 05:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')