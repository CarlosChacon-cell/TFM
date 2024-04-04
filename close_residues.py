


import pymol 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
args = parser.parse_args()


protein_name=args.protein.split('.')[0]
print(protein_name)

pymol.finish_launching()
cmd.load(args.protein)

cmd.select(f'sele, chain A w. 4 of chain B')
cmd.select('sele, br. sele')
residues='sele'

resset= set()
cmd.iterate(selector.process(residues), 'resset.add(resv)')

print(resset)

with open('residues.txt','w') as file:
    for element in resset:
        file.write(f'{int(element)-1}\n')

        
pymol.cmd.quit()






