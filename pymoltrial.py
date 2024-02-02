

#pymol distance getter. iterative process that checks if the distance between
#any aminoacid of the pirulo and any aminoacid of the interior of both proteins


# Import the PyMOL module
import pymol
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--peptide', '-pp', help='chain of the region we are changing', required=True)
parser.add_argument('--chains', '-ch' , help= 'chains that interact with the peptide', required=True)
parser.add_argument('--csv', '-csv', help='csv name, not .csv needed', required=True )
#parser.add_argument('--pept_res', '-pr', help='peptide residues we are interested in', required=True, nargs='+')
#parser.add_argument('--target_res', '-tr', help='target residues we are interested in', required=True, nargs= '+')
parser.add_argument('--distance_threshold', '-d', help='Distance Threshold. Default is 4', required=False, default=4)
parser.add_argument('--i', help='This is just teh iteration number because PMPNN_FR is very picky')
args = parser.parse_args()

# def generate_range(numbers):
#     start, end = map(int, numbers)
#     return list(range(start, end + 1))


# Initialize PyMOL
pymol.finish_launching()

# Load your structure
cmd.load(args.protein)

protein_name=args.protein.split('.')[0]
# List to store distances

dict={'pept_res':[ ], 'tar_res':[ ], 'distance':[ ], 'DA_AD': []}
interacting_residues=pd.DataFrame(dict)
# Residue indices
cmd.select(f'sele, chain {args.peptide}')
cmd.select('sele, br. sele')

peptide= "sele"
peptideres= set()
cmd.iterate(selector.process(peptide), 'peptideres.add(resv)')
peptide_residues = {args.peptide:peptideres}

#Trying to get the residue position from a selection

cmd.select(f'sele, chain {args.chains} w. 4 of chain {args.peptide}')
cmd.select('sele, br. sele')

interface = "sele"
storedresidues=set()
cmd.iterate(selector.process(interface), 'storedresidues.add(resv)')

#target indices 

target_residues = {args.chains:storedresidues}


# Calculate distance and store in the list.


for i in peptide_residues:
    for j in peptide_residues[i]:
        for k in target_residues: 
            for l in target_residues[k]:
                distance_polar = cmd.distance(f"chain {i} and resi {j} ", f"chain {k} and resi {l} ", mode=2)
                if distance_polar <4 and distance_polar != 0:
                    interacting_residues.loc[len(interacting_residues)]=[j,l,distance_polar, 'AD']
                    print ('residue added')
                else:
                    continue
                # distance_PP = cmd.distance(f"chain {i} and resi {j} ", f"chain {k} and resi {l} ", mode=5)
                # if distance_PP <4 and distance_PP != 0:
                #     interacting_residues.loc[len(interacting_residues)]=[j,l,distance_PP, 'AD']
                #     print ('residue added')
                # else:
                #     continue
print('Finished')
interacting_residues.to_csv(f'./{args.csv}.csv', sep='\t')

cmd.select(f'sele, {protein_name}')

cmd.save(f'output_{args.i}.pdb', 'sele')


pymol.cmd.quit()
