import numpy as np
import json

# Specify the path to your JSON file
file_path = '/emdata/cchacon/pae_trial.json'

# Open the JSON file in read mode
with open(file_path, 'r') as json_file:
    # Load the JSON data into a Python dictionary
    data = json.load(json_file)

# Now you can work with the 'data' dictionary containing your JSON content

pae=data['pae']

binderlen=439
filtered_pae=pae[binderlen:]
interacting_residues_pae=[]
residues=[]
mean=[]

with open('residues.txt', 'r') as resfile:
    for line in resfile:
        residues.append(int(line))


for lists in filtered_pae:
    for residue in residues:
        mean.append(lists[residue])
mean=np.mean(mean)