''''
This a code used for pairwise alignment between two protein sequences using the BlLOSUM62 matrix 
and a gap penalty of -10, as well as an extension penalty of -0.5.

Inputs:

--target: the original sequence 
--query: The sequence to be aligned 

Outputs:

pairwise_alignment.csv: A csv file with the number of identities between the sequences, the score of the alignment and the query and target name

'''



from Bio import Align
from Bio.Align import substitution_matrices 
from Bio import PDB
import pandas as pd
import argparse
import warnings
warnings.filterwarnings("ignore")

def extract_chain_sequence(pdb_file):
    # Create a parser object
    parser = PDB.PDBParser(QUIET=True)

    # Load the structure from the PDB file
    structure = parser.get_structure('protein', pdb_file)

    # Initialize an empty sequence string
    sequence = ''

    # Iterate over all models in the structure
    for model in structure:
        # Iterate over all chains in the model
        for chain in model:
            # Check if the chain id is 'A' (chain A)
            if chain.id == 'A':
                # Iterate over all residues in the chain
                for residue in chain:
                    # Check if the residue is an amino acid
                    if PDB.is_aa(residue):
                        # Get the residue's amino acid code and add it to the sequence
                        sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

    return sequence

def alignment(target_sequence, query_sequence, file_path='pairwise_alignment.csv'):

    target_sequence=extract_chain_sequence(args.target)
    query_sequence=extract_chain_sequence(args.query)

    aligner = Align.PairwiseAligner(open_gap_score = -10,extend_gap_score = -0.5 )

    matrix = substitution_matrices.load('BLOSUM62')

    aligner.substitution_matrix = matrix

    alignments= aligner.align(target_sequence, query_sequence)
    score = aligner.score (target_sequence, query_sequence)

    for alignment in alignments:
        c=alignment.counts()

    identities=c.identities

    dict= [{'score':score,
        'identities': identities,
        'query': args.query[:-4],
        'target': args.target[:-4]}]

    df=pd.DataFrame(dict)

    # Load existing CSV file into DataFrame
    try:
        existing_data = pd.read_csv(file_path)
        updated_data = pd.concat([existing_data,df], ignore_index=True)
        updated_data.to_csv(file_path, index=False)
    except FileNotFoundError:
        df.to_csv(file_path, index=False)

    print(f"Data appended and saved to {file_path}.")

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--target', help='target sequence')
    parser.add_argument('--query', help ='query sequence')
    args=parser.parse_args()
    alignment(args.target, args.query)


 