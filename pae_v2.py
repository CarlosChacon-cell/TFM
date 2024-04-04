import numpy as np
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--json', help='json file', required=True)
parser.add_argument('--protein', help='introduce your protein name, run_X_design_X')
args = parser.parse_args()

# Specify the path to your JSON file
file_path = args.json

#filepath=f'{args.protein}_dldesign0_pae.json'

# Open the JSON file in read mode
with open(file_path, 'r') as json_file:
    # Load the JSON data into a Python dictionary
    data = json.load(json_file)

# Now you can work with the 'data' dictionary containing your JSON content

pae=data['predicted_aligned_error']
residues=[]
mean=[]
residuefilename=f'{args.protein}_dldesign_0_cycle1_af2pred.txt'
with open('residues.txt', 'r') as resfile:
    for line in resfile:
        residues.append(int(line))

binderlen=residues[0]
residues=residues[1:]
filtered_pae=pae[binderlen:]

for lists in filtered_pae:
    for residue in residues:
        mean.append(lists[residue])
mean=np.mean(mean)

print('\n################\n')
print('\n################\n')
print('local_pae_interaction: ', mean)
print('\n################\n')
print('\n################\n')

pae1=pae[binderlen:]

pae2=pae[:binderlen]

pae1_filtered=[]
pae2_filtered=[]

for lista in pae1:
    pae1_filtered.append(lista[:binderlen])

for lista in pae2:
    pae2_filtered.append(lista[binderlen:])




pae_interaction=(np.mean(pae1_filtered)+np.mean(pae2_filtered))/2

print('\n################\n')
print('\n################\n')
print('global_pae_interaction: ', pae_interaction)
print('\n################\n')
print('\n################\n')
