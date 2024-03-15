
#This is just to build the ymla files neccesary for the RFAA from the PDB/fasta and compfile

import argparse
import re 

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help='protein fasta file path', required=False)
parser.add_argument('--pdb', help='protein pdb file path', required=False)
parser.add_argument('--comp', help='comp file path', required=True)
parser.add_argument('--folder', help='folder in which you want to save the yaml_file ', required=True)
args = parser.parse_args()


def pdb_to_fasta(pdb_file):
    sequence = ""

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[13:15] == "CA":
                # Extract the amino acid code
                aa = line[17:20].strip()

                # Append the amino acid to the sequence
                        

    return sequence




pattern=r'(?:.)(*$)'
comptype=re.search(pattern, args.comp)
comptype=comptype.group[0]
pattern1=r'(^*)(?:.)'
comp_name=re.search(pattern1, args.comp)
comp_name=comp_name.group[0]

if args.fasta:
    with open(f'{args.protein}', 'r') as fastafile:
        content = fastafile.read()
        # Split the content into lines
        lines = content.split('\n')
        # Iterate through the lines to find the protein name
        for line in lines:
            if line.startswith('>'):
                # Extract the protein name (assuming it's everything after '>')
                protein_name = line[1:].strip()
                break 
    filename=(f'{args.folder}/{protein_name}_{comp_name}')
    with open(filename, 'w') as file:
        file.write('defaults:\n'
                +'  '+'- base\n'
                + '\n'
                + '  ' + 'A:\n'
                + '\t'+f'fasta_file: {args.fasta}\n'
                + '\n'
                + 'sm_inputs:\n'
                + '  ' + 'B:\n'
                + '\t' + f'input: {args.com}'
                + '\t' + f'input_type: \"{comptype}\"')
                
if args.pdb:
    fasta_sequence=pdb_to_fasta(args.pdb)
    pattern = r'^[0-9a-zA-Z]'
    protein_name=re.search(pattern, args.pdb)
    protein_name=protein_name.group[0]

    filename=(f'{args.folder}/{protein_name}{comp_name}')

    with open(filename, 'w') as file:
        file.write('defaults:\n'
                +'  '+'- base\n'
                + '\n'
                + '  ' + 'A:\n'
                + '\t'+f'fasta_file: {args.fasta}\n'
                + '\n'
                + 'sm_inputs:\n'
                + '  ' + 'B:\n'
                + '\t' + f'input: {args.com}'
                + '\t' + f'input_type: \"{comptype}\"')


