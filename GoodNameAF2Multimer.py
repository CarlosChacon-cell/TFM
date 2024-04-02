
#No clue what for i wrote this code but seems unimportant
import pandas as pd
import os
import argparse

# Get the inputs
# parser = argparse.ArgumentParser()
# parser.add_argument('--target', '-t', help='target file path', required=True)
# parser.add_argument('--input', '-i', help='genes/proteins list file path', required=True)
# args = parser.parse_args()

# uniprot_df = pd.read_csv(args.input, na_values=['N/A']).dropna()
with open('./rcsb_pdb_7UXH', 'r') as target_fasta_file:
    fasta_lines = target_fasta_file.readlines()

# Create lists to store the data
fasta_names = []
fasta_sequences = []

# Parse the target.fasta data
for i in range(0, len(fasta_lines), 2):
    header = fasta_lines[i].strip().lstrip('>')
    sequence = fasta_lines[i + 1].strip()
    fasta_names.append(header)
    fasta_sequences.append(sequence)

# Create a dataframe from the parsed target.fasta data
fasta_df = pd.DataFrame({'FastaName': fasta_names, 'TargetSequence': fasta_sequences})

#Lets do the same with the ligands fasta

with open(args.input, 'r') as ligand_fasta_file:
    ligand_fasta_lines= ligand_fasta_file.readlines()

# Create lists to store the data
fasta_names = []
fasta_sequences = []

# Parse the target.fasta data
for i in range(0, len(ligand_fasta_lines), 2):
    header = ligand_fasta_lines[i].strip().lstrip('>')
    sequence = ligand_fasta_lines[i + 1].strip()
    fasta_names.append(header)
    fasta_sequences.append(sequence)

# Create a dataframe from the parsed target.fasta data
ligand_fasta_df = pd.DataFrame({'FastaName': fasta_names, 'TargetSequence': fasta_sequences})


