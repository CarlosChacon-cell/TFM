
#THis in reality is the instability index
import peptides as pep
import numpy as np
import re
import os

def pdb_to_fasta(pdb_file):
    sequence = ""

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[13:15] == "CA" and line[21] == "A":
                # Extract the amino acid code
                aa = line[17:20].strip()

                # Append the amino acid to the sequence
                sequence += aa

    return sequence


if __name__ == "__main__":
    # Get the input PDB file name from the user
    pdb_file = input("enter pdbfilepath")
    pattern=r'run_\d+'
    target_name= re.search(pattern, pdb_file)
    target_name=target_name.group(0)

    try:
        # Get the FASTA sequence
        fasta_sequence = pdb_to_fasta(pdb_file)
        peptide=pep.Peptide(fasta_sequence)
        instability_index=round(peptide.instability_index(), 2)
        hydrophobic_moment=round(np.mean(peptide.hydrophobic_moment_profile()), 3)
        print('\t'+ f'{instability_index}')
        print ('\t' +f'{hydrophobic_moment}')
        #We append the scores
        score_file="/emdata/cchacon/instability_index.csv"
        # Check if the file exists
        if not os.path.exists(score_file):
            # File doesn't exist, create it and write the header
             with open(score_file, 'w') as file:
                file.write('description'+ ','+ 'ins_index'+ ','+'hyd_moment' + '\n')

        with open(score_file, 'a') as file:
            file.write(f"{target_name},{instability_index},{hydrophobic_moment}\n")

    except Exception as e:
        print(f"Error: {e}")
