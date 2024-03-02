
#!/usr/bin/env python3

'''This code is just a small modification of CNIO code 
so the names of the FASTA output are legible and more easy to handle'''


import pandas as pd
import os
import argparse
import re


# Get the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--target', '-t', help='target file path', required=True)
parser.add_argument('--input', '-i', help='genes/proteins list file path', required=True)
args = parser.parse_args()

with open(args.target, 'r') as target_fasta_file:
    fasta_lines = target_fasta_file.readlines()

# Create lists to store the data
fasta_names = []
fasta_sequences = []
pattern1=r'^[^_]+'
# Parse the target.fasta data
for i in range(0, len(fasta_lines), 2):
    header = fasta_lines[i].strip().lstrip('>')
    target_name= re.search(pattern1, header)
    target_name=target_name.group(0)
    sequence = fasta_lines[i + 1].strip()
    fasta_names.append(target_name)
    fasta_sequences.append(sequence)

# Create a dataframe from the parsed target.fasta data
fasta_df = pd.DataFrame({'FastaName': fasta_names, 'TargetSequence': fasta_sequences})

#Lets do the same with the ligands fasta

with open( args.input,'r') as ligand_fasta_file:
    ligand_fasta_lines= ligand_fasta_file.readlines()

# Create lists to store the data
fasta_names = []
fasta_sequences = []

# Parse the target.fasta data
counter=1
for i in range(0, len(ligand_fasta_lines), 2):
    sequence = ligand_fasta_lines[i + 1].strip()
    if counter==1:
        fasta_names.append('OriginalSeq')
    else:
        fasta_names.append(f'NewSeq{counter-1}')
    fasta_sequences.append(sequence)
    counter=counter+1

# Create a dataframe from the parsed target.fasta data
ligand_fasta_df = pd.DataFrame({'LigandName': fasta_names, 'LigandSeq': fasta_sequences})

merged_df = pd.merge(fasta_df, ligand_fasta_df, how='cross')
#change the merged_df targetname column
merged_df['FastaName'] = merged_df['FastaName'] + '_' + merged_df['LigandName']
# Remove blank spaces in 'TargetSequence' and 'ProteinSequence'
merged_df['TargetSequence'] = merged_df['TargetSequence'].str.replace(' ', '')
merged_df['LigandSeq'] = merged_df['LigandSeq'].str.replace(' ', '')

# Create the new 'merged_sequences' column by combining 'TargetSequence' and 'ProteinSequence'
merged_df['merged_sequences'] = merged_df['TargetSequence'] + ':' + merged_df['LigandSeq']

# Create a new dataframe with 'FastaName' and 'merged_sequences' columns
new_df = merged_df[['FastaName', 'merged_sequences']].copy()
# OUTPUT
output_directory = 'AF2_runs'
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
pattern1=r'^[^_]+'

for index, row in new_df.iterrows():
    fasta_name = row['FastaName']
    merged_sequences = row['merged_sequences']
    fasta_folder = os.path.join(output_directory, fasta_name, 'input')
    os.makedirs(fasta_folder, exist_ok=True)
    target_fasta_file_path = os.path.join(fasta_folder, f'{fasta_name}.fasta')
    with open(target_fasta_file_path, 'w') as target_fasta_file:
        target_fasta_file.write(f'>{fasta_name}\n{merged_sequences}')


print("Folders and input files created successfully.")

# AF2 Submission commands
slurm_script = "af2_submit_script"
output_script_filename = "submit_af2_screen.sh"
with open(output_script_filename, "w") as output_script_file:
    for index, row in new_df.iterrows():
        fasta_name = row['FastaName']
        job_command = f"sbatch -t 1400 {slurm_script} {os.path.abspath(os.path.join(output_directory, fasta_name))}\n"
        output_script_file.write(job_command)

print(f"SLURM submission script saved to {output_script_filename}")



