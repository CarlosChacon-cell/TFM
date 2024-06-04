


import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

df=pd.read_csv("/Users/carlo/Downloads/HOMOLUMO.csv", index_col=0)

df['IP']=-df['HOMO']*1.3023+0.481
df['EA']=-df['LUMO']*0.6091-0.475

metrics_df=pd.read_csv("/Users/carlo/Downloads/LigandSearchMetrics.csv", index_col=0,header=0, names=['Compound', 'mean_plddt', 'pae_interaction'])

merge_df=pd.merge(df, metrics_df, how='inner', on='Compound')

# plt.figure(figsize=(10,8))
# sns.scatterplot(data=merge_df, x='Gap', y='IP', hue='pae_interaction', s=100, palette='crest')
# plt.show()


# plt.figure(figsize=(12,10))
# sns.set_theme(style='darkgrid')
# sns.scatterplot(data=merge_df, x='pae_interaction', y='mean_plddt', s=100)
# plt.xlabel(f'PAE Interaction ($\AA$)')
# plt.ylabel(f'pLDDT')
# plt.title('pLDDT vs PAE interaction')
# plt.show()

df_energy=pd.read_csv("/Users/carlo/Downloads/ETOT.csv", header=0)

merge_df=pd.merge(merge_df,df_energy, on='Compound', how='inner')
merge_df.to_csv('/Users/carlo/Desktop/quantic.csv')

# plt.figure(figsize=(8,10))
# sns.set_theme(style='darkgrid')
# sns.scatterplot(data=merge_df, x='ETOT', y='Gap', palette='crest',s=100)
# plt.xlabel('ETOT (eV)')
# plt.ylabel('HOMO-LUMO gap (eV)')
# plt.title('ETOT vs HOMO-LUMO gap')
# plt.show()

# fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(8,10))
# sns.set_theme(style='darkgrid')
# ax1.set_title('ETOT Distribution')
# fig_1=sns.violinplot(data=merge_df, y='ETOT',ax=ax1)
# fig_1.invert_yaxis()
# ax1.set_ylabel('ETOT (eV)')
# ax2.set_title('HLG Distribution')
# sns.violinplot(data=merge_df, y='Gap',ax=ax2)
# ax2.set_ylabel('HLG (eV)')
# plt.show()

# plt.figure(figsize=(10,8))
# sns.set_theme(style='darkgrid')
# sns.scatterplot(data=merge_df, x='EA', y='IP', hue='Gap', s=100, palette='crest')
# plt.xlabel('Electron Affinity (eV)')
# plt.ylabel('Ionization Potential (eV)')
# plt.legend(title='HOMO-LUMO Gap')
# plt.title('EA vs IP For Predicted Ligands')
# plt.show()

plt.figure(figsize=(10,8))
sns.set_theme(style='darkgrid')
sns.scatterplot(data=merge_df, x='IP', y='pae_interaction', hue='Gap', s=100, palette='crest')
plt.xlabel('Ionization Potential (eV)')
plt.ylabel(f'PAE interaction ($\AA$)')
plt.legend(title='HOMO-LUMO Gap')
plt.title('IP vs PAE For Predicted Ligands')
plt.show()