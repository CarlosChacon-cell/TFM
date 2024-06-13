import sys
from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help = "input file")
parser.add_argument("--chain", "-c", help="chain to extract")
args, unknown = parser.parse_known_args()

def remove_remark_lines(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if not line.startswith('REMARK'):
                f_out.write(line)

def extract_chain_sequence(input_file, chain_id):
    sequences = {}
    with open(input_file, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            chain = record.annotations['chain']
            if chain == chain_id:
                sequences[chain] = record.seq
    return sequences

pattern = r'run_\d+_design_\d+.*_dldesign_(\d+).pdb'
PDBFile = args.input
PDB_out = 'noremarks.pdb'

try:
    fasta_name = re.search(pattern, PDBFile).group(0)
except AttributeError:
    fasta_name = PDBFile[:-4]
remove_remark_lines(PDBFile, PDB_out)
if args.chain:
    chain=args.chain
    chain_sequences = extract_chain_sequence(PDB_out, chain)

fasta_file = f'{fasta_name}_chain{args.chain}.fasta'
with open(fasta_file, 'w') as fasta:
    for chain_id, seq in chain_sequences.items():
        fasta.write(f'>{fasta_name}_{chain_id}\n{seq}\n')

print(f"Chain A sequence extracted and saved to '{fasta_file}'")
