import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re 

# Define a function to calculate the Tanimoto similarity between two SMILES strings
def calculate_tanimoto(smiles1, smiles2):
    # Convert SMILES strings to RDKit Mol objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    # Generate molecular fingerprints (e.g., Morgan fingerprints, also known as ECFP)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=4096)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=4096)
    
    # Calculate the Tanimoto similarity
    tanimoto_similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    return tanimoto_similarity

def calculate_similarity_matrix(similarity):
    n = len(similarity)
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim_matrix[i, j] = similarity[i]
            sim_matrix[j, i] = similarity[i]
    return 1 - sim_matrix #Distance Matrix

def clustering(distance_matrix):
    # Perform Agglomerative Clustering
    cluster = AgglomerativeClustering(n_clusters=100, affinity='precomputed', linkage='complete')
    cluster.fit(distance_matrix)

    # Output cluster labels
    labels = cluster.labels_
    print('Cluster labels:', labels)
    sns.clustermap(distance_matrix, row_cluster=True, col_cluster=True, figsize=(8, 8))
    plt.show()
    return cluster

if __name__ == '__main__' :
    pattern=r'>\s+<smiles>*'
    dict_tc={
        'Comp1':[],
        'Comp2':[],
        'Tc':[]
    }
    interaction_df=pd.read_csv('LigandSearchMetrics.csv',header=0)
    filtered_df=interaction_df[interaction_df['pae_interaction'] < 10]
    comp_name=[]
    for comp in filtered_df['Comp_name']:
        comp_name.append(comp)
    for h in range(len(comp_name)-1):
        comp1=comp_name[h]
        with open(f'/emdata/cchacon/RFAA/20240520_throughput_myb/compounds/{comp1}.sdf', 'r') as target_file:
            lines=target_file.readlines()
            i=0
            for line in lines:
                if re.match(pattern,line):
                    target_smile=lines[i+1]
                    break
                else:
                    i+=1
        for j in range(h+1,len(comp_name)):
            comp2=comp_name[j]
            with open(f'/emdata/cchacon/RFAA/20240520_throughput_myb/compounds/{comp2}.sdf', 'r') as query_file:
                lines=query_file.readlines()
                i=0
                for line in lines:
                    if re.match(pattern,line):
                        query_smile=lines[i+1]
                        break
                    else:
                        i+=1
            similarity = calculate_tanimoto(target_smile, query_smile)
            dict_tc['Comp1'].append(comp1)
            dict_tc['Comp2'].append(comp2)
            dict_tc['Tc'].append(round(similarity,2))
    df=pd.DataFrame(dict_tc)
    df.to_csv('Tanimoto_similarity_hits.csv')
    distance_matrix=calculate_similarity_matrix(dict_tc['Tc'])
    cluster=clustering(distance_matrix)
    print(cluster)


# Convert similarity matrix to distance matrix (1 - similarity)

