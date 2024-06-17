'''
This code computes which residues are in the interacting interface (Distance between proteins < 4 A) between ligand and protein.
--> Input:
    --protein: protein file path, in pdb format
    --length: Length of the protein target. This is needed since the pdb output from RFAA numbers each atom from the small molecule as a residue
               of the target.
-->Output:
    close_residues.csv: CSV file with the protein name, its length, the interacting surface in %, the interacting residues of the target 
'''




import argparse
import os
import pymol

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--length', '-l', help='The length of the protein', type=int)
args = parser.parse_args()
protein_name=args.protein.split('/')[0]
try:
    protein_name=protein_name.split('.')[0]
except:
    protein_name=protein_name




pymol.finish_launching()
cmd.load(args.protein)

# Select interacting residues
cmd.delete('sele')
cmd.select(f'sele, resi 1-{args.length} w. 4 of resn LG1')
cmd.select('sele, br. sele')
residues_sel = 'sele'

res_set = set()
cmd.iterate(residues_sel, 'res_set.add(resv)')

#compute interacting surface

interacting_surface=len(res_set)/args.length*100

# Save data to CSV
filename = f'close_residues.csv'
file_exists = os.path.isfile(filename)

if file_exists:
    with open(filename, 'a') as file:
        # Convert the list of residues to a whitespace-separated string representation
        res_str = ' '.join(map(str, list(res_set)))
        file.write(f'{protein_name},{args.length},{interacting_surface},{res_str}\n')
else:
    with open(filename, 'w') as file:
        # Write column headers
        file.write('compound_name,length,interacting_surface,interacting_residues\n')
        # Convert the list of residues to a whitespace-separated string representation
        res_str = ' '.join(map(str, list(res_set)))
        file.write(f'{protein_name},{args.length},{interacting_surface},{res_str},\n')

pymol.cmd.quit()
