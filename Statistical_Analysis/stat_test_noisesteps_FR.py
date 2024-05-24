from scipy.stats import chi2_contingency
import pandas as pd 
import numpy as np
from scipy.stats import ttest_ind
import argparse


df=pd.read_csv('/emdata/cchacon/RFD_PD_test/output_af2.csv')
# Define the data
positives_FR = 14
total_FR = 755
positives_no_FR = 10
total_no_FR = 1290

# Create a contingency table
observed = [[positives_FR, total_FR - positives_FR], [positives_no_FR, total_no_FR - positives_no_FR]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic:", chi2)
print("P-value:", p)


# Data
pae_interaction_FR = df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='campaign_run_112_FR')]
pae_interaction_no_FR = df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='campaign_run_112_noFR')]
plddt_binder_FR = df['pae_interaction'][(df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_FR')]
plddt_binder_no_FR = df['pae_interaction'][(df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_noFR')]

# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_FR, pae_interaction_no_FR, equal_var=False)
print(f'T-Test for pae_interaction between FR and no FR: t-stat = {t_stat}, p-value = {p_value}')

# T-Test for plddt_binder
t_stat, p_value = ttest_ind(plddt_binder_FR, plddt_binder_no_FR, equal_var=False)
print(f'T-Test for plddt_binder: t-stat = {t_stat}, p-value = {p_value}')


pae_interaction_ns5=df['pae_interaction'][(df['pae_interaction'] < 15) & (df['campaign']=='campaign_run_112_ns5')]
pae_interaction_ns10=df['pae_interaction'][(df['pae_interaction'] < 15) & (df['campaign']=='campaign_run_112_ns10')]
pae_interaction_ns20=df['pae_interaction'][(df['pae_interaction'] < 15) & (df['campaign']=='campaign_run_112_ns20')]
pae_interaction_ns30=df['pae_interaction'][(df['pae_interaction'] < 15) & (df['campaign']=='campaign_run_112_ns30')]

bonferroni_correction=0.05/6
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns5, pae_interaction_ns10, equal_var=False)
print(f'T-Test for pae_interaction between 5 and 10 noise steps: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns10, pae_interaction_ns20, equal_var=False)
print(f'T-Test for pae_interaction between 10 and 20 noise steps: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns20, pae_interaction_ns30, equal_var=False)
print(f'T-Test for pae_interaction between 20 and 30 noise steps: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns20, pae_interaction_ns5, equal_var=False)
print(f'T-Test for pae_interaction between 20 and 5 noise steps: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns5, pae_interaction_ns30, equal_var=False)
print(f'T-Test for pae_interaction between 30 and 5 noise steps: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns10, pae_interaction_ns30, equal_var=False)
print(f'T-Test for pae_interaction between 10 and 30 noise steps: t-stat = {t_stat}, p-value = {p_value}')

positives_5ns=len(df['pae_interaction'][(df['pae_interaction'] < 10) & (df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_ns5')])
total_5ns=len(df['pae_interaction'][(df['campaign']=='campaign_run_112_ns5')])

positives_10ns=len(df['pae_interaction'][(df['pae_interaction'] < 10) & (df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_ns10')])
total_10ns=len(df['pae_interaction'][(df['campaign']=='campaign_run_112_ns10')])

positives_20ns=len(df['pae_interaction'][(df['pae_interaction'] < 10) & (df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_ns20')])
total_20ns=len(df['pae_interaction'][(df['campaign']=='campaign_run_112_ns20')])

positives_30ns=len(df['pae_interaction'][(df['pae_interaction'] < 10) & (df['plddt_binder'] > 80) & (df['campaign']=='campaign_run_112_ns30')])
total_30ns=len(df['pae_interaction'][(df['campaign']=='campaign_run_112_ns30')])


# Create a contingency table
observed = [[positives_5ns, total_5ns - positives_5ns], [positives_10ns, total_10ns- positives_10ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 5 and 10 ns:", chi2)
print("P-value between 5 and 10 ns:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')

# Create a contingency table
observed = [[positives_20ns, total_20ns - positives_20ns], [positives_5ns, total_5ns- positives_5ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 5 and 20:", chi2)
print("P-value between 5 and 20:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_5ns, total_5ns- positives_5ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 5 and 30:", chi2)
print("P-value between 5 and 30:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_10ns, total_10ns- positives_10ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 10 and 30:", chi2)
print("P-value between 10 and 30:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')

# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_20ns, total_20ns- positives_20ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 20 and 30:", chi2)
print("P-value between 20 and 30:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
# Create a contingency table
observed = [[positives_20ns, total_20ns - positives_20ns], [positives_10ns, total_10ns- positives_10ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 10 and 20:", chi2)
print("P-value between 10 and 20:", p)
if p < bonferroni_correction:
    print('########\n')
    print('This difference is significant\n')
    print('#########\n')
