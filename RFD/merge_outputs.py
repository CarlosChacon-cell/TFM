

import pandas as pd
import os
import glob 


def pdb_exists(pdb_id, folder):
    pdb_file = os.path.join(folder, f"{pdb_id}.pdb")
    return os.path.exists(pdb_file)


# Define paths to the CSV files
folder1 ='/emdata/cchacon/RFD_partial_diff/'
folder2 = '/emdata/cchacon/RFD_PD_test'
filename = 'output_af2.csv'

# Load CSV files
csv1 = pd.read_csv(os.path.join(folder1, filename))
csv2 = pd.read_csv(os.path.join(folder2, filename))

pdbs=glob.glob('/emdata/cchacon/RFD_partial_diff/*pdb')
filtered_csv1=pd.DataFrame(columns=csv1.columns)
for pdb in pdbs:
    pdb=pdb.split('/')[4]
    pdb=pdb.split('.')[0]
    filtered_csv1=pd.concat([filtered_csv1,(csv1.loc[(csv1['description']==pdb) & (csv1['pae_interaction']<10)])])

pdbs2=glob.glob('/emdata/cchacon/RFD_PD_test/*pdb')
filtered_csv2 = pd.DataFrame(columns=csv2.columns)
for pdb in pdbs2:
    pdb=pdb.split('/')[4]
    pdb=pdb.split('.')[0]
    filtered_csv2=pd.concat([filtered_csv2,(csv2.loc[(csv2['description']==pdb) & (csv2['pae_interaction']<10)])])


# Merge CSV files
merged_csv = pd.concat([filtered_csv1, filtered_csv2], ignore_index=True)

# Save merged CSV file
output_folder = '/emdata/cchacon/'
output_filename = 'merged_output_pd.csv'
merged_csv.to_csv(os.path.join(output_folder, output_filename), index=False)

print("Merge complete. Merged file saved at:", os.path.join(output_folder, output_filename))
