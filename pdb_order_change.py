#!/usr/bin/env python3


'''This code changes the order of the PDB, putting in first place the chain we are interested in
Probabily this has been stolen out of CNIO code but seems written by chatGPT so maybe 
I'm the "author" of this '''


import argparse
from Bio.PDB import PDBParser, PDBIO

parser=argparse.ArgumentParser()
parser.add_argument('--input', '-i' , required=True, help='input file')
parser.add_argument('--output', '-o', required=True)
parser.add_argument('--chain', '-c', required=True, type=str)
args=parser.parse_args()





# Replace 'input.pdb' and 'output.pdb' with the actual input and output PDB file paths
input_pdb_file = args.input
output_pdb_file = args.output

# Specify the target chain ID you want to move to the beginning
target_chain_id = args.chain

# Use Biopython's PDBParser to parse the input PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', input_pdb_file)

# Get a list of chain IDs in the original order
original_chain_order = [chain.id for chain in structure[0]]

# Move the target chain to the beginning of the chain order
new_chain_order = [target_chain_id] + [chain_id for chain_id in original_chain_order if chain_id != target_chain_id]

# Create a new structure with the modified chain order
new_structure = parser.get_structure('modified_protein', input_pdb_file)
new_structure[0].child_list = [structure[0][chain_id] for chain_id in new_chain_order]

# Use Biopython's PDBIO to write the modified structure to a new PDB file
io = PDBIO()
io.set_structure(new_structure)
io.save(output_pdb_file)
