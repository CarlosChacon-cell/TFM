#!/usr/bin/env python3

import numpy as np
import json
import argparse
import pandas as pd
import re 
import glob

parser = argparse.ArgumentParser()
parser.add_argument('--protein', type=str, help='introduce your protein name, run_X_design_X')
args = parser.parse_args()

# Specify the path to your JSON file
pattern_run = r'run_\d+'
pattern_design = r'design_\d+'

protein_run_match = re.search(pattern_run, args.protein)
protein_design_match = re.search(pattern_design, args.protein)

protein_run = protein_run_match.group() if protein_run_match else None
protein_design = protein_design_match.group() if protein_design_match else None


try:
    file_path = glob.glob(f'../output/{protein_run}/pae*{protein_run}*{protein_design}*.json')[0]
    protein_name=args.protein.split('.')[0]
except IndexError:
    try:
        protein_name=args.protein.split('.')[0]
        file_path=glob.glob(f'../output/{protein_run}/pae_{protein_run}_{protein_design}*json')[0]
    except:
        protein_name=args.protein.split('.')[0]
        file_path=glob.glob(f'../jsons/pae_{protein_run}_{protein_design}*json')[0]



# Open the JSON file in read mode
with open(file_path, 'r') as json_file:
    # Load the JSON data into a Python dictionary
    data = json.load(json_file)

# Now you can work with the 'data' dictionary containing your JSON content

pae=data['predicted_aligned_error']
plddt=data['plddt']
residues=[]
mean=[]
residuefilename='close_residues.csv'

residues_df=pd.read_csv(residuefilename, index_col=False)

binderlen=residues_df['length'][residues_df['protein_name']==protein_name]
binderlen=int(binderlen.iloc[0])
interacting_surface=float(residues_df['interacting_surface'][residues_df['protein_name']==protein_name].iloc[0])
binder_residues=residues_df['binder_residues'][residues_df['protein_name']==protein_name].str.split()
binder_residues=binder_residues.to_list()[0]
target_residues=residues_df['target_residues'][residues_df['protein_name']==protein_name].str.split()
target_residues=target_residues.to_list()[0]


if type(residues)==float:
    pae1=pae[binderlen:]

    pae2=pae[:binderlen]

    pae1_filtered=[]
    pae2_filtered=[]

    for lista in pae1:
        pae1_filtered.append(lista[:binderlen])

    for lista in pae2:
        pae2_filtered.append(lista[binderlen:])

    pae_interaction_global=(np.mean(pae1_filtered)+np.mean(pae2_filtered))
    pae_interaction_local=pae_interaction_global
else:
    binder_residues=[int(number)-1 for number in binder_residues]
    target_residues=[int(number)-1 for number in target_residues]
    filtered_pae=pae[binderlen:]

    for residue_target in target_residues:
        lists=pae[residue_target]
        for binder_residue in binder_residues:
            mean.append(lists[binder_residue]/(plddt[residue_target]/100))
    for residue_binder in binder_residues:
        lists=pae[residue_binder]
        for residue_target in target_residues:
            mean.append(lists[residue_target]/(plddt[residue_binder]/100))

    pae_interaction_local=np.mean(mean)
    if pae_interaction_local > 30.0:
        pae_interaction_local=30.0

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


file_path='pae_local_global.csv'
# Load existing CSV file into DataFrame
try:
    existing_data = pd.read_csv('pae_local_global.csv')
    updated_data = pd.concat([existing_data,df], ignore_index=True)
    updated_data.to_csv('pae_local_global.csv', index=False)
except FileNotFoundError:
    df.to_csv(file_path, index=False)

print(f"Data appended and saved to {file_path}.")