
# This is the new section to compute posible pipi interactions. This is not a code in itself but just a trial to check if the ideas worked. I have to implement this into pymol trial, but i will do that when cation-pi also works
import pymol
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--peptide', '-pp', help='chain of the region we are changing', required=True)
parser.add_argument('--chains', '-ch' , help= 'chains that interact with the peptide', required=True)
parser.add_argument('--csv', '-csv', help='csv name, not .csv needed', required=True )
#parser.add_argument('--pept_res', '-pr', help='peptide residues we are interested in', required=True, nargs='+')
#parser.add_argument('--target_res', '-tr', help='target residues we are interested in', required=True, nargs= '+')
parser.add_argument('--distance_threshold', '-d', help='Distance Threshold. Default is 4', required=False, default=4)
# parser.add_argument('--i', help='This is just teh iteration number because PMPNN_FR is very picky')
args = parser.parse_args()

# def generate_range(numbers):
#     start, end = map(int, numbers)
#     return list(range(start, end + 1))

#compute the angle between the rings using the dot product
def vecAngle(vec1, vec2):
    '''

    return the angle between them in degree.
    '''
    dotprod = vec1[0][0]*vec2[0][0]+vec1[0][1]*vec2[0][1]+vec1[0][2]*vec2[0][2]
    magnitude1=np.sqrt(vec1[0][0]**2+vec1[0][1]**2+vec1[0][2]**2)
    magnitude2=np.sqrt(vec2[0][0]**2+vec2[0][1]**2+vec2[0][2]**2)
    deg = np.arccos(round(dotprod/(magnitude1*magnitude2), 4)) * 180 / np.pi
    print(f'Angle between centroids is {deg}')
    if deg > 90:
        deg = 180 - deg
    return deg


# Initialize PyMOL
pymol.finish_launching()

# Load your structure
cmd.load(args.protein)

protein_name=args.protein.split('.')[0]
# List to store distances

dict={'pept_res':[ ], 'tar_res':[ ], 'distance':[ ], 'Polar_PP': []}
interacting_residues=pd.DataFrame(dict)
                
# This first block search for ligand residues that can be involved in the interaction, restricting the search to only aromatic residues
cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and (resn PHE or resn TYR or resn TRP or resn HIS) w. 20 of chain {args.chains} and (resn PHE or resn TYR or resn TRP or resn HIS)')
cmd.select('sele, br. sele')
ligand_rings = "sele"
storedligrings=set()
cmd.iterate(selector.process(ligand_rings), 'storedligrings.add(resv)') #This is how pymol works, no clue about why

lig_aro_residues={args.peptide:storedligrings}
print(lig_aro_residues)
#Just the same with the protein chains 
cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn PHE or resn TYR or resn TRP) w. 20 of chain {args.peptide} and (resn PHE or resn TYR or resn TRP)')
cmd.select('sele, br. sele')
protein_rings = "sele"
storedprotrings=set()
cmd.iterate(selector.process(protein_rings), 'storedprotrings.add(resv)')

prot_aro_residues={args.chains:storedprotrings}
print(prot_aro_residues)

#we define the thresholds

dih_parallel=25
dih_tshape=80

#A loop to generate the pseudoatoms at the centroids of the residues 

for i in lig_aro_residues:
    for j in lig_aro_residues[i]:
        lig_center=cmd.pseudoatom(f'lig_center_{j}', f'{protein_name} and chain {args.peptide} and resi {j} and name cg+cz')
for k in prot_aro_residues: 
            for l in prot_aro_residues[k]:
                prot_center=cmd.pseudoatom(f'prot_center_{l}', f'{protein_name} and chain {args.chains} and resi {l} and name cg+cz')


# Another loop to store this residues 
for i in lig_aro_residues:
    for j in lig_aro_residues[i]:
        for k in prot_aro_residues: 
            for l in prot_aro_residues[k]:
                distance_pp = cmd.distance(f'lig_center_{j}////ps1', f'prot_center_{l}////ps1')
                print(f'distance is {distance_pp} \n')
                xyz_lig=cmd.get_coords(f'lig_center_{j}////ps1')
                print(f'{xyz_lig} \n')
                xyz_prot=cmd.get_coords(f'prot_center_{l}////ps1')
                angle=vecAngle(xyz_lig, xyz_prot)
                if (distance_pp <15 and distance_pp != 0.0) and (angle < dih_parallel or angle > dih_tshape) :
                    interacting_residues.loc[len(interacting_residues)]=[j,l,distance_pp, 'PP']
                    print ('residue added')


interacting_residues.to_csv(f'./{args.csv}.csv', sep='\t')