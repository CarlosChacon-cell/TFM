
#Open Pymol, select chain A and save it on the new folder 

import pymol
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--folder')
args = parser.parse_args()

pattern=r'run_\d+'
target_name= re.search(pattern, args.protein)
target_name=target_name.group(0)


pymol.finish_launching()

cmd.load(args.protein)
cmd.select("ligand, chain A")


cmd.save(f'{args.folder}/{target_name}.pdb', 'sele')
pymol.cmd.quit()