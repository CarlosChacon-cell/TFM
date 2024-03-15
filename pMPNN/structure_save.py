import pymol 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--i', help='This is just teh iteration number because PMPNN_FR is very picky')
parser.add_argument('--folder', help='Folder to save the file')
args = parser.parse_args()

# Initialize PyMOL
pymol.finish_launching()

# Load your structure
cmd.load(args.protein)

protein_name=args.protein.split('.')[0]

#Save your structure
cmd.select(f'sele, {protein_name}')
if args.folder:
    cmd.save(f'{args.folder}/protein_{args.i}.pdb', 'sele')
else: 
    cmd.save(f'protein_{args.i}.pdb', 'sele')


print(f'output_{args.i}.pdb')



pymol.cmd.quit()