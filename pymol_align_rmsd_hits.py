#!/usr/bin/env python3

import pymol
from pymol import cmd
import argparse
import glob

### Parsing input file info
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help = "input folder")
parser.add_argument("--og", "-o", help="Original input")
args, unknown = parser.parse_known_args()

if not args.input:
    print("Please provide an input folder to process")
    exit()
else:
    input = str(args.input)

print(f'Running pymol on {input} ...')

# Load the reference and target protein structures
prefix=input

pdb_list=glob.glob(f'*.pdb') #Meant to be run from a hits folder or something similar
protein_names=[]


for pdb in pdb_list:
    pdb=pdb.split('.')[0]
    protein_names.append(pdb)
i=0 
print(protein_names)   
for pdb in pdb_list:
    cmd.load(pdb,f"{protein_names[i]}")
    i+=1

cmd.load(args.og,"original")

for protein in protein_names:
    cmd.align(f"{protein}","original")


#Calculate RMSD between chains
chains = ["A","B"]
csv_file_name = f'{input}/rmsd.csv'
with open(csv_file_name, 'w') as csv_file:
    csv_file.write("Rank,Chain,Global RMSD\n")
    for protein in protein_names:
        global_rmsd=cmd.rms_cur(f'{protein} and chain A and name CA' , 'original and chain A and name CA')
        #selective_rmsd = cmd.rms_cur(f'{protein} and chain A and b>50', f'original and chain A and b>50')
        csv_file.write(f"{protein},'A',{global_rmsd}\n")

print(f'RMSD values have been saved to {csv_file_name}')

# for protein in protein_names:
#     cmd.save(f'{protein}_aln.pdb', f'{protein}')


print("Done!")