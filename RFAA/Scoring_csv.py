import argparse
import pandas as pd
from Scoring import extract_scoring, extract_features
import glob 
import numpy as np

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder', '-f')
    args=parser.parse_args()
    csv_name='close_residues.csv'
    filename=glob.glob(f'{args.folder}/*_aux.pt')[0]
    mean_plddt, pae_interaction, pae, plddt=extract_scoring(filename)
    #PAE is a pytorch.Tensor with only two dimensions, so we remove the tensor structure to make things easier
    pae=pae[0]
    plddt=plddt[0]
    #Extract the several things
    protein_length, interacting_surface, residues= extract_features(csv_name, filename)
    compound_length=len(pae)- protein_length 
    #We extract the cutre pae into one list
    cutre=[]
    for residue in residues:
        for atom in range(protein_length,len(pae)):
            cutre.append(float(pae[int(residue)][atom]/plddt[int(residue)]))
            cutre.append(float(pae[atom][int(residue)]/plddt[atom]))
    cutre_local=np.mean(cutre)
    score_dict={
        'Comp_name':args.folder,
        'mean_plddt': [round(mean_plddt,3)],
        'pae_interaction': [round(pae_interaction,3)],
        'pae_interaction_local':[round(cutre_local,3)]
    }
    score_df=pd.DataFrame(score_dict)
    try:
        old_df=pd.read_csv('LigandSearchMetrics.csv', index_col=0)
        merge_df=pd.concat([old_df, score_df])
        merge_df.to_csv('LigandSearchMetrics.csv')
    except:
        score_df.to_csv('LigandSearchMetrics.csv')
        
