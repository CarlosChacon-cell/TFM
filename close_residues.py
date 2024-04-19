import argparse
import os
import pymol

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
args = parser.parse_args()

protein_name = args.protein.split('.')[0]

pymol.finish_launching()
cmd.load(args.protein)

# First compute the length
cmd.select('sele, chain A')
cmd.select('sele, br. sele')
all_sel = 'sele'
all_set = set()
cmd.iterate(all_sel, 'all_set.add(resv)')

length = len(all_set)

# Select interacting residues
cmd.delete('sele')
cmd.select('sele, chain A w. 4 of chain B')
cmd.select('sele, br. sele')
residues_sel = 'sele'

res_set = set()
cmd.iterate(residues_sel, 'res_set.add(resv)')

#compute interacting surface

interacting_surface=len(res_set)/length*100

# Save data to CSV
filename = f'close_residues.csv'
file_exists = os.path.isfile(filename)

if file_exists:
    with open(filename, 'a') as file:
        # Convert the list of residues to a whitespace-separated string representation
        res_str = ' '.join(map(str, list(res_set)))
        file.write(f'{protein_name},{length},{interacting_surface},{res_str}\n')
else:
    with open(filename, 'w') as file:
        # Write column headers
        file.write('protein_name,length,interacting_surface,interacting_residues\n')
        # Convert the list of residues to a whitespace-separated string representation
        res_str = ' '.join(map(str, list(res_set)))
        file.write(f'{protein_name},{length},{interacting_surface},{res_str},\n')

pymol.cmd.quit()
