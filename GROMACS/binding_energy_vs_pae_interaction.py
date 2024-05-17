import csv
import pandas as pd
import argparse
import numpy as np
import os
import glob
import re


def extract_protein_name(pdb_input):
    pattern=r'input/(.*)\.pdb'
    protein_name=re.search(pattern, pdb_input).group(1)
    return protein_name

# Define a function to extract values from the final column
def extract_total_values(input_file, output_file):
    start_processing = False  # Flag to indicate when to start processing rows
    with open(input_file, 'r') as csv_input, open(output_file, 'w', newline='') as csv_output:
        reader = csv.reader(csv_input)
        writer = csv.writer(csv_output)
        for row in reader:
            # Check if the row contains the "Delta Energy Terms" string
            if 'Delta Energy Terms' in ','.join(row):
                start_processing = True
                continue  # Skip this row and move to the next
            # Start processing rows below the "Delta Energy Terms" line
            if start_processing:
                # Check if the row has the correct number of columns
                if len(row) >= 13:
                    # Extract and store the value from the final column
                    total_value = row[-1]
                    writer.writerow([total_value])

# Call the function to extract and store total values

def extract_pae_interaction(input_pdb):
    with open(input_pdb, 'r') as file:
        for line in file:
            if line.startswith('pae_interaction'):
                pae_interaction=line[-6:-1]
    return pae_interaction


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder', help='Folder in which the pdb and the delta are contained')
    args=parser.parse_args()
    try:
        input_file=glob.glob(os.path.join(args.folder, '*csv'))[0]
        output_file=os.path.join(args.folder, 'delta.csv')
        input_pdb=glob.glob(os.path.join('input/*pdb'))[0]
        extract_total_values(input_file, output_file)
        df=pd.read_csv(output_file, index_col=None, header=0)
        print(df)

        delta=np.mean(df)
        pae_interaction=extract_pae_interaction(input_pdb)
        protein_name=extract_protein_name(input_pdb)
        print(protein_name)
        with open('../deltavspae.csv', 'a') as file:
            file.write(f'{protein_name},{delta},{pae_interaction}\n')

    except:
        print('Nope')