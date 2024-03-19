
#This is just to build the ymla files neccesary for the RFAA from the PDB/fasta and compfile
import os 
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
                sequence+=aa

                # Append the amino acid to the sequence
                        
    return sequence


def get_file_type(file_path):
    _, file_extension = os.path.splitext(file_path)
    file_extension=file_extension.split('.')[1]
    return file_extension.lower()

def get_file_name(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]





comp_path = args.comp
comptype = get_file_type(comp_path)
comp_name = get_file_name(comp_path)


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
                + '\t' + f'input: {args.comp}\n'
                + '\t' + f'input_type: \"{comptype}\"')
                
if args.pdb:
    fasta_sequence=pdb_to_fasta(args.pdb)
    protein_name=get_file_name(args.pdb)
    fasta_name = f'fasta/{protein_name}.fasta'
    with open(fasta_name, 'w') as fasta:
        fasta.write(f'>{protein_name}\n' + 
                   fasta_sequence)
        
    filename=(f'{args.folder}/{protein_name}_{comp_name}.yaml')

    with open(filename, 'w') as file:
        file.write('defaults:\n'
                +'  '+'- base\n'
                + '\n'
                + '  ' + 'A:\n'
                + '\t'+f'fasta_file: {fasta_name}\n'
                + '\n'
                + 'sm_inputs:\n'
                + '  ' + 'B:\n'
                + '\t' + f'input: {args.comp}\n'
                + '\t' + f'input_type: \"{comptype}\"')


