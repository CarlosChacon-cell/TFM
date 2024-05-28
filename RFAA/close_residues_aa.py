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
