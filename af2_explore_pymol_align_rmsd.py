#!/usr/bin/env python3

import pymol
from pymol import cmd
import argparse
import glob

### Parsing input file info
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help = "input folder", required=True)
args, unknown = parser.parse_known_args()

if not args.input:
    print("Please provide an input folder to process")
    exit()
else:
    input = str(args.input)

print(f'Running pymol on {input} ...')

# Load the reference and target protein structures
prefix=input

pdb_list=glob.glob('*.pdb') #Meant to be run from a hits folder or something similar

for pdb in pdb_list:
    cmd.load(pdbin1[0], "rank1")


cmd.align("rank2 and chain A and b>50","rank1 and chain A and b>50")
cmd.align("rank3 and chain A and b>50","rank1 and chain A and b>50")
cmd.align("rank4 and chain A and b>50","rank1 and chain A and b>50")
cmd.align("rank5 and chain A and b>50","rank1 and chain A and b>50")

# Calculate RMSD between chains
ranks = ["rank2", "rank3", "rank4", "rank5"]
chains = ["A", "B", "C"]
csv_file_name = f'{input}/{prefix}_rmsd.csv'
with open(csv_file_name, 'w') as csv_file:
    csv_file.write("Rank,Chain,Global RMSD,Selective RMSD (B>50)\n")
    for rank in ranks:
        for chain in chains:
            global_rmsd = cmd.rms_cur(f"{rank} and chain {chain}", f"rank1 and chain {chain}")
            selective_rmsd = cmd.rms_cur(f"{rank} and chain {chain} and b>50", f"rank1 and chain {chain} and b>50")
            csv_file.write(f"{rank},{chain},{global_rmsd},{selective_rmsd}\n")

print(f'RMSD values have been saved to {csv_file_name}')

# Save the aligned models
pdbout1=f'{pdbin1[0][:-4]}_aln.pdb'
pdbout2=f'{pdbin2[0][:-4]}_aln.pdb'
pdbout3=f'{pdbin3[0][:-4]}_aln.pdb'
pdbout4=f'{pdbin4[0][:-4]}_aln.pdb'
pdbout5=f'{pdbin5[0][:-4]}_aln.pdb'

cmd.save(pdbout1, "rank1")
cmd.save(pdbout2, "rank2")
cmd.save(pdbout3, "rank3")
cmd.save(pdbout4, "rank4")
cmd.save(pdbout5, "rank5")

print("Done!")