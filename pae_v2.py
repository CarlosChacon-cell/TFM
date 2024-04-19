#!/usr/bin/env python3

import numpy as np
import json
import argparse
import pandas as pd 
parser = argparse.ArgumentParser()
parser.add_argument('--protein', help='introduce your protein name, run_X_design_X')
args = parser.parse_args()

# Specify the path to your JSON file
file_path = f'./jsons/pae_{args.protein[5:-4]}.json'
protein_name=args.protein.split('.')[0]


# Open the JSON file in read mode
with open(file_path, 'r') as json_file:
    # Load the JSON data into a Python dictionary
    data = json.load(json_file)

# Now you can work with the 'data' dictionary containing your JSON content

pae=data['predicted_aligned_error']
residues=[]
mean=[]
residuefilename='/emdata/cchacon/TRF1/20240415_campaign_hotspotDIM_TRF1/close_residues.csv'

residues_df=pd.read_csv(residuefilename, index_col=False)

binderlen=residues_df['length'][residues_df['protein_name']==protein_name]
binderlen=int(binderlen.iloc[0])
interacting_surface=float(residues_df['interacting_surface'][residues_df['protein_name']==protein_name].iloc[0])
residues=residues_df['interacting_residues'][residues_df['protein_name']==protein_name].str.split()
residues=residues.to_list()[0]
if type(residues)==float:
    pae_interaction_global=35
    pae_interaction_local=35
else:
    residues=[int(number)-1 for number in residues]
    filtered_pae=pae[binderlen:]

    for lists in filtered_pae:
        for residue in residues:
            mean.append(lists[residue])
    pae_interaction_local=np.mean(mean)

    pae1=pae[binderlen:]

    pae2=pae[:binderlen]

    pae1_filtered=[]
    pae2_filtered=[]

    for lista in pae1:
        pae1_filtered.append(lista[:binderlen])

    for lista in pae2:
        pae2_filtered.append(lista[binderlen:])

    pae_interaction_global=(np.mean(pae1_filtered)+np.mean(pae2_filtered))/2


df=pd.DataFrame({
'pae_interaction_local':[pae_interaction_local], 
'pae_interaction_global':[pae_interaction_global],
'protein_name':[args.protein.split('.')[0]],
'length': [binderlen],
'interacting_surface':[interacting_surface]
})

print('\n################\n')
print('\n################\n')
print('local_pae_interaction: ', pae_interaction_local)
print('\n################\n')
print('\n################\n')

print('\n################\n')
print('\n################\n')
print('global_pae_interaction: ', pae_interaction_global)
print('\n################\n')
print('\n################\n')
file_path='pae_local_global.csv'
# Load existing CSV file into DataFrame
try:
    existing_data = pd.read_csv('pae_local_global.csv')
    updated_data = pd.concat([existing_data,df], ignore_index=True)
    updated_data.to_csv('pae_local_global.csv', index=False)
except FileNotFoundError:
    df.to_csv(file_path, index=False)

print(f"Data appended and saved to {file_path}.")