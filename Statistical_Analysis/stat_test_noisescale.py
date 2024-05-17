from scipy.stats import chi2_contingency
import pandas as pd 
import numpy as np
from scipy.stats import ttest_ind
import argparse


df=pd.read_csv('/emdata/cchacon/partial_diff_noise_scale_test/output_af2.csv')

pae_interaction_ns005=df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='run_112_noise_scale_005')]
pae_interaction_ns01=df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='run_112_noise_scale_01')]
pae_interaction_ns05=df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='run_112_noise_scale_05')]
pae_interaction_ns1=df['pae_interaction'][(df['pae_interaction'] < 10) & (df['campaign']=='run_112_noise_scale_1')]


# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns01, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 01 noise scale: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 05 noise scale: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns005, pae_interaction_ns1, equal_var=False)
print(f'T-Test for pae_interaction between 005 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns01, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 01 and 05 noise scale: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns01, pae_interaction_ns1, equal_var=False)
print(f'T-Test for pae_interaction between 01 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')
# T-Test for pae_interaction
t_stat, p_value = ttest_ind(pae_interaction_ns1, pae_interaction_ns05, equal_var=False)
print(f'T-Test for pae_interaction between 05 and 1 noise scale: t-stat = {t_stat}, p-value = {p_value}')

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

# Create a contingency table
observed = [[positives_20ns, total_20ns - positives_20ns], [positives_5ns, total_5ns- positives_5ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 5 and 20:", chi2)
print("P-value between 5 and 20:", p)

# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_5ns, total_5ns- positives_5ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 5 and 30:", chi2)
print("P-value between 5 and 30:", p)

# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_10ns, total_10ns- positives_10ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 10 and 30:", chi2)
print("P-value between 10 and 30:", p)


# Create a contingency table
observed = [[positives_30ns, total_30ns - positives_30ns], [positives_20ns, total_20ns- positives_20ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 20 and 30:", chi2)
print("P-value between 20 and 30:", p)

# Create a contingency table
observed = [[positives_20ns, total_20ns - positives_20ns], [positives_10ns, total_10ns- positives_10ns]]

# Perform chi-square test
chi2, p, _, _ = chi2_contingency(observed)

print("Chi-square statistic between 10 and 20:", chi2)
print("P-value between 10 and 20:", p)
