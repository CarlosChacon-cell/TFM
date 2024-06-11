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

# Select binder interacting residues
cmd.delete('sele')
cmd.select('sele, chain A w. 4 of chain B')
cmd.select('sele, br. sele')
binder_sel = 'sele'


binder_set = set()
cmd.iterate(binder_sel, 'binder_set.add(resv)')

#compute binder interacting surface

interacting_surface=len(binder_set)/length*100

# Select target interacting residues
cmd.delete('sele')
cmd.select('sele, chain B w. 4 of chain A')
cmd.select('sele, br. sele')
target_sel = 'sele'


target_set = set()
cmd.iterate(target_sel, 'target_set.add(resv)')

# Save data to CSV
filename = f'close_residues.csv'
file_exists = os.path.isfile(filename)

if file_exists:
    with open(filename, 'a') as file:
        # Convert the list of residues to a whitespace-separated string representation
        binder_str = ' '.join(map(str, list(binder_set)))
        target_str=' '.join(map(str, list(target_set)))
        file.write(f'{protein_name},{length},{interacting_surface},{binder_str},{target_str}\n')
else:
    with open(filename, 'w') as file:
        # Write column headers
        file.write('protein_name,length,interacting_surface,binder_residues,target_residues\n')
        # Convert the list of residues to a whitespace-separated string representation
        binder_str = ' '.join(map(str, list(binder_set)))
        target_str=' '.join(map(str, list(target_set)))
        file.write(f'{protein_name},{length},{interacting_surface},{binder_str},{target_str}\n')
pymol.cmd.quit()
