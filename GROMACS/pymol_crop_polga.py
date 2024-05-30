import pymol 
import glob 
import re

protein_files=glob.glob('*modified*.pdb')
pattern=r'(PolGa_mut\d+).*'
for protein in protein_files:
    print(protein)
    protein_name=re.search(pattern, protein).group(1)
    pymol.finish_launching()
    cmd.load(protein)
    cmd.delete('loop')
    cmd.select('loop',f'resi 991-1052')
    cmd.save (f'{protein_name}_loop.pdb','loop')
    cmd.delete('all')
pymol.cmd.quit()


