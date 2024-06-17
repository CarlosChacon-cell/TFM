#!/usr/bin/env python3

'''
This is a code to compute the RMSD between the original input of the partial diffusion and the hits generated. 
It considers the RMSD of the C alpha.

Inputs:

--input: Folder in which there is the PDBs (I recommend running this code from the folder itself)
--og: The original pdb fom which we are running the partial diffusion 

Output:

--rmsd.csv: A csv file with the RMSD, the query and the original name and the chain we are using to compute the RMSD (in our case always A so it is hardcoded)

WARNING! This code must be run inside PyMol: pymol -c /path/to/pymol_align_rmsd_hits.py --input example.pdb --og original.pdb --chain A
'''




import pymol
from pymol import cmd
import argparse
import glob

### Parsing input file info
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help = "input folder", default='.')
parser.add_argument("--og", "-o", help="Original input")
parser.add_argument('--chain', '-c', help='Chain to compute the RMSD', default='A', type=str)
args, unknown = parser.parse_known_args()

if not args.input:
    print("Please provide an input folder to process")
    exit()
else:
    input = str(args.input)

pymol.finish_launching()



# Load the reference and target protein structures
prefix=input

pdb_list=glob.glob(f'*.pdb') #Meant to be run from a hits folder or something similar
protein_names=[]


for pdb in pdb_list:
    pdb=pdb.split('.')[0]
    protein_names.append(pdb)
i=0 
for pdb in pdb_list:
    cmd.load(pdb,f"{protein_names[i]}")
    i+=1

cmd.load(args.og,"original")

for protein in protein_names:
    cmd.align(f"{protein} and chain B","original and chain B")


#Calculate RMSD between chains
chains = ["A","B"]
csv_file_name = f'{input}/rmsd.csv'
with open(csv_file_name, 'w') as csv_file:
    csv_file.write("query,Chain,Global RMSD\n")
    for protein in protein_names:
        global_rmsd=cmd.rms_cur(f'{protein} and chain A and name CA' , f'original and chain A and name CA')
        #selective_rmsd = cmd.rms_cur(f'{protein} and chain A and b>50', f'original and chain A and b>50')
        csv_file.write(f"{protein},'A',{global_rmsd}\n")

print(f'RMSD values have been saved to {csv_file_name}')

# for protein in protein_names:
#     cmd.save(f'{protein}_aln.pdb', f'{protein}')


print("Done!")

pymol.cmd.quit()
