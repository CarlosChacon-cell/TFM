

import numpy as np
import os
import argparse
import re
import pandas as pd

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--pdbdir", type=str)
parser.add_argument("--pdbs", type=str)
parser.add_argument("--indices", required=False, nargs='+')
parser.add_argument("--csv", required = False)
parser.add_argument("--verbose", action="store_true", default=False)
args = parser.parse_args()

#A function to generate a range of position from 2 numbers 
def generate_range(numbers):
    start, end = map(int, numbers)
    return list(range(start, end + 1))

#Now i put the indices, generating a range if there is any hyphen
indices=[]
if args.indices:
    positions=args.indices
    positions=positions[0].split(',')
    pattern= r'\d+-\d+'
    for number in positions:
        match = re.search(pattern, number)
        if match:
            numbers=list(number.split('-'))
            numbers=generate_range(numbers)
            for number in numbers:
                indices.append(number)
        else:
            indices.append(int(number))

#read the interacting.csv and add new residues to the indices
if args.csv:
    interacting_df=pd.read_csv(f'./{args.csv}.csv', sep='\t')
    for i in interacting_df['pept_res']:
        if i not in indices:
            indices.append(int(i))
        else:
            continue


#pdb path
pdb_path=args.pdbs

print(f"Adding FIXED labels to {pdb_path} at positions {indices}")

remarks = []
for position in indices:
    remark = f"REMARK PDBinfo-LABEL:{position: >5} FIXED"
    remarks.append(remark)

# Uncomment below to add hotspots
remarks_str = '\n'.join(remarks)
with open(pdb_path, 'a') as f:
    f.write('\n')
    f.write(remarks_str)