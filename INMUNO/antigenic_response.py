'''
Based on this strange paper https://febs.onlinelibrary.wiley.com/doi/epdf/10.1016/0014-5793%2890%2980535-Q
They get an aminoacid surface propensity and an aminoacid antigen propensity
then they calculate the average propensity of an heptamere as :

<A>=1/N * sum(Ai) , where Ai is fsi/fagi . This vlaue is assigned
to the middle antigen in the 7mer. They do this for the whole petide
Then, they compute the <Aprotein>. If it is larger than
1, any 7mer with residues above 1 are antigenic (or bigger stretch)
If it is smaller than 1, any 7mer residues above the avg
are antigenic

This is only valid to compute the probability of a linear epitope for B-cells, so its usefulness is very limited
'''



import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import glob
import argparse




def antigenic_calculator(file_path):
    #Input variables
    dict_antigen={
    'A':1.064,
    'C':1.412,
    'D':0.866,
    'E':0.851,
    'F':1.091,
    'G':0.874,
    'H':1.105,
    'I':1.152,
    'K':0.930,
    'L':1.250,
    'M':0.826,
    'N':0.776,
    'P':1.064,
    'Q':1.015,
    'R':0.873,
    'S':1.012,
    'T':0.909,
    'V':1.383,
    'W':0.893,
    'Y':1.161
    }

    with open (file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                protein_name=line.strip('>')
                continue
            else:
                sequence=line

    #output variables
    out_dict= {
        'protein_name':[],
        'Antigenic_propensity':[],
        'N_antigenic_seq':[]
    }


    antigenic_values=[]
    for i in range(len(sequence)-7):
        value=[]
        window=sequence[i:i+7]
        for letter in window:
            value.append(dict_antigen[letter])
        antigenic_values.append(np.mean(value))

    protein_mean_antigenic=np.mean(antigenic_values)
    print('The mean antigenic value of your protein is: ', protein_mean_antigenic)


    i=0
    while i < (len(antigenic_values)-8):
        if protein_mean_antigenic > 1:
            if all(j > 1 for j in antigenic_values[i:i+8]):
                print('Possible antigenic sequence: ', sequence[i+3:i+11], ' at position ', i+3)
                i+=8
            i+=1
        else:
            if all(j > protein_mean_antigenic for j in antigenic_values[i:i+8]):
                print('Possible antigenic sequence: ', sequence[i+3:i+11])
                i+=8
            i+=1 

    plt.plot(antigenic_values, label='Antigenic value')
    plt.title('Antigenic Value vs Residue Position')
    plt.xlabel('Antigenic value')
    plt.ylabel('Residue number-3')
    plt.hlines(y=protein_mean_antigenic, xmin=0, xmax=len(sequence)-3, color='r', label = 'Mean antigenic value')
    plt.legend()
    # plt.show()
    return protein_name, protein_mean_antigenic


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder', '-f', help='Folder in which the fasta are stored')
    args=parser.parse_args()
    antigenic_dict={
        'protein_name':[],
        'mean_antigenic':[]
    }
    file_list= glob.glob(f'{args.folder}/*fasta')
    for file_path in file_list:
        protein_name, mean_antigenic = antigenic_calculator(file_path)
        antigenic_dict['protein_name'].append(protein_name.strip())
        antigenic_dict['mean_antigenic'].append(mean_antigenic)
    df=pd.DataFrame(antigenic_dict)
    df.to_csv('antigenic_values.csv', index=False)
    

