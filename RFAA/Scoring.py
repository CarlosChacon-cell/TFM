'''
Code that provides the needed function for the CUTRE_AA.py to work
'''
import torch 
import argparse
import pandas as pd 
import re 
import numpy as np

def extract_features(csv_name, pt_name):
    pattern=r'(ZINC\d+)_aux.pt'
    compound_name=re.search(pattern, pt_name).group(1)
    residues_df=pd.read_csv(csv_name, index_col=False, header=0)
    protein_lenght=int(residues_df['length'][residues_df['compound_name']==compound_name].iloc[0])
    interacting_surface=float(residues_df['interacting_surface'][residues_df['compound_name']==compound_name].iloc[0])
    residues=residues_df['interacting_residues'][residues_df['compound_name']==compound_name].str.split()
    residues=residues.to_list()[0]
    return protein_lenght, interacting_surface, residues

def extract_scoring(filename):
    scoring_dict=torch.load(filename, map_location='cpu')
    mean_plddt=scoring_dict['mean_plddt']
    pae_interaction=scoring_dict['pae_inter']
    pae=scoring_dict['pae']
    plddt=scoring_dict['plddts']
    return mean_plddt, pae_interaction, pae, plddt


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='File to be analysed')
    args=parser.parse_args()
    filename=args.file
    csv_name='close_residues.csv'
    mean_plddt,pae_interaction,pae, plddt=extract_scoring(filename)
    #PAE is a pytorch.Tensor with only two dimensions, so we remove the tensor structure to make things easier
    pae=pae[0]
    #Extract the several things
    protein_length, interacting_surface, residues= extract_features(csv_name, filename)
    compound_length=len(pae)- protein_length 
    plddt=plddt[0]
    #We extract the cutre pae into one list
    cutre=[]
    for residue in residues:
        for atom in range(protein_length,len(pae)):
            cutre.append(float(pae[int(residue)][atom]/plddt[int(residue)]))
            cutre.append(float(pae[atom][int(residue)]/plddt[atom]))
    print(cutre)
    cutre_local=np.mean(cutre)
    print('The mean plddt of the prediction is: ', mean_plddt)
    print('\nThe global pae_interaction is: ', pae_interaction)
    print('\n The CUTRE score is: ', cutre_local)