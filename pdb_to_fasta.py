import sys
from Bio import SeqIO
import re

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

pattern = r'run_\d+'
PDBFile = sys.argv[1]
PDB_out = f'{PDBFile[:-4]}_noremarks.pdb'
fasta_name = re.search(pattern, PDBFile).group(0)

remove_remark_lines(PDBFile, PDB_out)
chain_sequences = extract_chain_sequence(PDB_out, 'A')

fasta_file = f'{fasta_name}_chainA.fasta'
with open(fasta_file, 'w') as fasta:
    for chain_id, seq in chain_sequences.items():
        fasta.write(f'>{fasta_name}_{chain_id}\n{seq}\n')

print(f"Chain A sequence extracted and saved to '{fasta_file}'")
