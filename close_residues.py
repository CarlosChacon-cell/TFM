
import pymol 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
args = parser.parse_args()


protein_name=args.protein.split('.')[0]
print(protein_name)

pymol.finish_launching()
cmd.load(args.protein)

#First compute the length 
cmd.select('sele, chain A')
cmd.select('sele, br. sele')
all='sele'
allset= set()
cmd.iterate(selector.process(all), 'allset.add(resv)')

length=len(allset)

filename=f'{protein_name}_residues.txt'

with open(filename,'w') as file:
    file.write(f'{length}\n')
cmd.delete('sele')
cmd.select(f'sele, chain A w. 4 of chain B')
cmd.select('sele, br. sele')
residues='sele'

resset= set()
cmd.iterate(selector.process(residues), 'resset.add(resv)')

print(resset)

with open(filename,'a') as file:
    for element in resset:
        file.write(f'{int(element)}\n')

pymol.cmd.quit()






