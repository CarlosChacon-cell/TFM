
# This is the new section to compute posible pi-cation interactions. 
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
cmd.select(f'sele, chain {args.peptide} and (resn PHE or resn TYR or resn TRP) w. 10 of chain {args.chains} and (resn ARG or resn LYS)  and name N* ')
cmd.select('sele, br. sele')
ligand_rings = "sele"
storedligrings=set()
cmd.iterate(selector.process(ligand_rings), 'storedligrings.add(resv)') #This is how pymol works, no clue about why

lig_aro_residues={args.peptide:storedligrings}
print(lig_aro_residues)

#Just the same with the protein chains 
cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn PHE or resn TYR or resn TRP) w. 10 of chain {args.peptide} and (resn ARG or resn LYS)  and name N*')
cmd.select('sele, br. sele')
protein_rings = "sele"
storedprotrings=set()
cmd.iterate(selector.process(protein_rings), 'storedprotrings.add(resv)')

prot_aro_residues={args.chains:storedprotrings}


# This second block search for ligand positive charges that can be involved in the cation pi interaction. If we follow the bibliography, this search should be restricted to the amino groups of R and K. But every N group could be involved ?
cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and (resn ARG or resn LYS)  and name NH* w. 10 of chain {args.chains} and (resn PHE or resn TYR or resn TRP)')
cmd.select('sele, br. sele')
ligand_cations = "sele"
storedligcations=set()
cmd.iterate(selector.process(ligand_cations), 'storedligcations.add(resv)') #This is how pymol works, no clue about why

lig_cat_residues={args.peptide:storedligcations}

#Just the same with the protein chains 
cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn ARG or resn LYS)  and name NH* w. 10 of chain {args.peptide} and (resn PHE or resn TYR or resn TRP)')
cmd.select('sele, br. sele')
protein_cations = "sele"
storedprotcations=set()
cmd.iterate(selector.process(protein_cations), 'storedprotcations.add(resv)')

prot_cat_residues={args.chains:storedprotcations}

#we define the thresholds

dih_parallel=25
dih_tshape=80

### TRP are problematic, so maybe i have to do a little trick to compute its centroid. As i have read, the centroid must be located on the six C ring, and this centroids seems ot be located between CE3 and CZ2. So the pseudoatom of the TRP must be computed differently.
###What im gonna do is to do a new selection with the resi of the TRP and then if the number of the other selections coincides with this one I calculate teh pseudoatom ina different manner.

#Selecting just the Trps of ligand

cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and resn TRP ')
cmd.select('sele, br. sele')
lig_trp = "sele"
storedligtrp=set()
cmd.iterate(selector.process(lig_trp), 'storedligtrp.add(resv)')

lig_trp_residues={args.peptide:storedligtrp}

#selecting the Trps of the protein

cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and resn TRP ')
cmd.select('sele, br. sele')
prot_trp = "sele"
storedprottrp=set()
cmd.iterate(selector.process(prot_trp), 'storedprottrp.add(resv)')

prot_trp_residues={args.chains:storedprottrp}

#A loop to generate the pseudoatoms at the centroids of the residues 

for i in lig_aro_residues:
    for j in lig_aro_residues[i]:
        if j not in lig_trp_residues[args.peptide]:
            lig_center=cmd.pseudoatom(f'lig_center_{j}', f'{protein_name} and chain {args.peptide} and resi {j} and name cg+cz')
        else:
             lig_center=cmd.pseudoatom(f'lig_center_{j}', f'{protein_name} and chain {args.peptide} and resi {j} and name ce3+cz2')

for k in prot_aro_residues: 
            for l in prot_aro_residues[k]:
                if l not in prot_trp_residues[args.chains]:
                    prot_center=cmd.pseudoatom(f'prot_center_{l}', f'{protein_name} and chain {args.chains} and resi {l} and name cg+cz')
                else:
                    prot_center=cmd.pseudoatom(f'prot_center_{l}', f'{protein_name} and chain {args.chains} and resi {l} and name ce3+cz2')


# Another loop to store this residues 
for i in lig_aro_residues:
    for j in lig_aro_residues[i]:
        for k in prot_cat_residues[args.chains]:
                distance_cp = cmd.distance(f'lig_center_{j}////ps1', f'chain {args.chains} and resi {k} and name NH*')
                print(f'distance is {distance_cp} \n')
               #Angle im not sure how to compute it or which angles to use 
                if (distance_cp < 8 and distance_cp != 0.0) :
                    interacting_residues.loc[len(interacting_residues)]=[j,l,distance_cp, 'CP']
                    print ('residue added')
for k in prot_aro_residues: 
    for l in prot_aro_residues[k]:
                for j in lig_cat_residues[args.peptide]:
                    distance_cp = min(cmd.distance(f' chain {args.peptide} and resi {j} and name NH*', f'prot_center_{l}////ps1'))
                    print(f'distance is {distance_cp} \n')
                    if (distance_cp < 8 and distance_cp != 0.0):
                        interacting_residues.loc[len(interacting_residues)]=[j,l,distance_cp, 'CP']
                        print ('residue added')


interacting_residues.to_csv(f'./{args.csv}.csv', sep='\t')