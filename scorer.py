
#V1.2
#Removing protein_name beacause it is not needed and fuck the process
#20240305
#pymol distance getter. iterative process that checks if the distance between
#any aminoacid of the pirulo and any aminoacid of the interior of both proteins


# Import the PyMOL module
import pymol
import argparse
import numpy as np
import os


parser = argparse.ArgumentParser()
parser.add_argument('--protein', '-pt', help='protein file path', required=True)
parser.add_argument('--peptide', '-pp', help='chain of the region we are changing', required=True)
parser.add_argument('--chains', '-ch' , help= 'chains that interact with the peptide', required=True)
parser.add_argument('--csv', '-csv', help='csv name, not .csv needed')
#parser.add_argument('--pept_res', '-pr', help='peptide residues we are interested in', required=True, nargs='+')
#parser.add_argument('--target_res', '-tr', help='target residues we are interested in', required=True, nargs= '+')
parser.add_argument('--distance_threshold', '-d', help='Distance Threshold. Default is 4', required=False, default=4)
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
    print(dotprod)
    magnitude1=np.sqrt(vec1[0][0]**2+vec1[0][1]**2+vec1[0][2]**2)
    print(magnitude1)
    magnitude2=np.sqrt(vec2[0][0]**2+vec2[0][1]**2+vec2[0][2]**2)
    deg = np.arccos(round(dotprod/(magnitude1*magnitude2), 4)) * 180 / np.pi
    print(deg)
    if deg > 90:
        deg = 180 - deg
    return deg

def pairwise_dist(sel1, sel2, max_dist, output="N", sidechain="N", show="N"):
    cmd.delete("dist*")
    extra=""
    if sidechain=="Y":
        extra=" and not name c+o+n"

    #builds models
    m1 = cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
    m1o = cmd.get_object_list(sel1)
    m2 = cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
    m2o = cmd.get_object_list(sel2)

    #defines selections
    cmd.select("__tsel1a", sel1+" around "+str(max_dist)+" and "+sel2+extra)
    cmd.select("__tsel1", "__tsel1a and "+sel2+extra)
    cmd.select("__tsel2a", sel2+" around "+str(max_dist)+" and "+sel1+extra)
    cmd.select("__tsel2", "__tsel2a and "+sel1+extra)
    cmd.select("IntAtoms_"+str(max_dist), "__tsel1 or __tsel2")
    cmd.select("IntRes_"+str(max_dist), "byres IntAtoms_"+str(max_dist))

    #controlers-1
    if len(m1o)==0: 
        print("warning, '"+sel1+extra+"' does not contain any atoms.")
        return
    if len(m2o)==0: 
        print("warning, '"+sel2+extra+"' does not contain any atoms.")
        return
    distances=[]
    #measures distances
    s=""
    counter=0
    for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
            distances.append(math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord,m2.atom[c2].coord)))))
            distance=min(distances)
            if distance<float(max_dist):
                s+="%s/%s/%s/%s/%s to %s/%s/%s/%s/%s: %.3f\n" % (m1o[0],m1.atom[c1].chain,m1.atom[c1].resn,m1.atom[c1].resi,m1.atom[c1].name,m2o[0],m2.atom[c2].chain,m2.atom[c2].resn,m2.atom[c2].resi,m2.atom[c2].name, distance)
                counter+=1
                if show=="Y": cmd.distance (m1o[0]+" and "+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name, m2o[0]+" and "+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)

    #controler-2
    if counter==0: 
        return

    #outputs
    if output=='A':
        return distances 
    if output=="S":
        print(s)
    if output=="P":
        f=open('IntAtoms_'+max_dist+'.txt','w')
        f.write("Number of distances calculated: %s\n" % (counter))
        f.write(s)
        f.close()
        print("Results saved in IntAtoms_%s.txt" % max_dist)
    print("Number of distances calculated: %s" % (counter))
    cmd.hide("lines", "IntRes_*")
    if show=="Y": cmd.show("lines","IntRes_"+max_dist)
    cmd.deselect()

cmd.extend("pairwise_dist", pairwise_dist)









#global variable
protein_name=args.protein.split('.')[0]
print(protein_name)

dict={'pept_res':[ ], 'tar_res':[ ], 'distance':[ ], 'type': [], 'score':[]}

# Initialize PyMOL
pymol.finish_launching()

# Load your structure
cmd.load(args.protein)

protein_name=args.protein.split('.')[0]
# List to store distances

# Residue indices
cmd.select(f'sele, chain {args.peptide} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his) w. 4 of chain {args.chains} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his)')
peptide='sele'
peptideres= set()
cmd.iterate(selector.process(peptide), 'peptideres.add(resv)')
peptide_residues = {args.peptide:peptideres}
print(peptide_residues)
#Trying to get the residue position from a selection
cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his) w. 4 of chain {args.peptide} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his)')
interface = "sele"
storedresidues=set()
cmd.iterate(selector.process(interface), 'storedresidues.add(resv)')
    
#target indices 

target_residues = {args.chains:storedresidues}
print(target_residues)

#Calculate distance and store in the list.


for i in peptide_residues:
    for j in peptide_residues[i]:
        for k in target_residues: 
            for l in target_residues[k]:
                distance_polar=pairwise_dist(f'chain {i} and resi {j}',f'chain {k} and resi {l}', max_dist=4, output='A')
                if distance_polar != None:
                    distance_polar=min(distance_polar)
                    dict['pept_res'].append(j)
                    dict["tar_res"].append(l)
                    dict['distance'].append(distance_polar)
                    dict['type'].append("Polar")
                    dict['score'].append(-0.25)
                    print ('residue added')
                else:
                     continue


#add the selector for hydrophobes in teh peptide

cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and (resn ala+gly+val+ile+leu+phe+met) w. 5 of chain {args.chains} and (resn ala+gly+val+ile+leu+phe+met)')
cmd.select(f'sele, br. sele')
hydro_ligand="sele"
storedhydro_ligand=set()
cmd.iterate(selector.process(hydro_ligand), 'storedhydro_ligand.add(resv)')

ligand_hydro_residues={args.peptide:storedhydro_ligand}



#Select the hydrophobes of the target

cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn ala+val+ile+leu+phe+met) w. 5 of chain {args.peptide} and (resn ala+val+ile+leu+phe+met)')
cmd.select(f'sele, br. sele')
hydro_target="sele"
storedhydro_target=set()
cmd.iterate(selector.process(hydro_target), 'storedhydro_target.add(resv)')

target_hydro_residues={args.chains:storedhydro_target}

#time for the calculus

for i in ligand_hydro_residues:
    for j in ligand_hydro_residues[i]:
        for k in target_hydro_residues: 
            for l in target_hydro_residues[k]:
                distance_hydro=pairwise_dist(f'chain {i} and resi {j}',f'chain {k} and resi {l}', max_dist=4, output='A')
                if distance_hydro != None:
                    distance_hydro=min(distance_hydro)
                    print('residue added')
                    dict['pept_res'].append(j)
                    dict["tar_res"].append(l)
                    dict['distance'].append(distance_hydro)
                    dict['type'].append("H-H")
                    dict['score'].append(-1)
                else:
                    continue

#Now, hydrophobic ligand, polar target interaction
                
cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and (resn ala+gly+val+ile+leu+phe+met) w. 5 of chain {args.chains} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his)')
cmd.select(f'sele, br. sele')
hydro_polar_ligand="sele"
storedhydro_polar_ligand=set()
cmd.iterate(selector.process(hydro_polar_ligand), 'storedhydro_polar_ligand.add(resv)')

ligand_hp_residues={args.peptide:storedhydro_polar_ligand}

#Same with the other peptides 

cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his) w. 5 of chain {args.peptide} and (resn ala+gly+val+ile+leu+phe+met)')
cmd.select(f'sele, br. sele')
polar_hydro_target="sele"
storedpolar_hydro_target=set()
cmd.iterate(selector.process(polar_hydro_target), 'storedpolar_hydro_target.add(resv)')

target_ph_residues={args.chains:storedpolar_hydro_target}

#time for the calculus

for i in ligand_hp_residues:
    for j in ligand_hp_residues[i]:
        for k in target_ph_residues: 
            for l in target_ph_residues[k]:
                distance_hp=pairwise_dist(f'chain {i} and resi {j}',f'chain {k} and resi {l}', max_dist=4, output='A')
                if distance_hp != None:
                    distance_hp=min(distance_hp)
                    print('residue added')
                    dict['pept_res'].append(j)
                    dict["tar_res"].append(l)
                    dict['distance'].append(distance_hp)
                    dict['type'].append("H-P")
                    dict['score'].append(0.6)
                else:
                    continue


#Now, hydrophobic target, polar ligand interaction
                
cmd.delete('sele')
cmd.select(f'sele, chain {args.chains} and (resn ala+gly+val+ile+leu+phe+met) w. 5 of chain {args.peptide} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his)')
cmd.select(f'sele, br. sele')
hydro_polar_target="sele"
storedhydro_polar_target=set()
cmd.iterate(selector.process(hydro_polar_target), 'storedhydro_polar_target.add(resv)')

target_hp_residues={args.chains:storedhydro_polar_target}

#Same with the other peptides 

cmd.delete('sele')
cmd.select(f'sele, chain {args.peptide} and (resn ser+thr+cys+tyr+asn+gln+asp+glu+lys+arg+his) w. 5 of chain {args.chains} and (resn ala+gly+val+ile+leu+phe+met)')
cmd.select(f'sele, br. sele')
polar_hydro_ligand="sele"
storedpolar_hydro_ligand=set()
cmd.iterate(selector.process(polar_hydro_ligand), 'storedpolar_hydro_ligand.add(resv)')

ligand_ph_residues={args.peptide:storedpolar_hydro_ligand}

#time for the calculus

for i in ligand_ph_residues:
    for j in ligand_ph_residues[i]:
        for k in target_hp_residues: 
            for l in target_hp_residues[k]:
                distance_ph=pairwise_dist(f'chain {i} and resi {j}',f'chain {k} and resi {l}', max_dist=4, output='A')
                if distance_ph != None:
                    distance_ph=min(distance_ph)
                    dict['pept_res'].append(j)
                    dict["tar_res"].append(l)
                    dict['distance'].append(distance_ph)
                    dict['type'].append("P-H")
                    dict['score'].append(0.6)
                else:
                    continue

print('Finished')

output_file = f'{args.csv}.csv'
total_score=0
for i in dict['score']:
     print(total_score)
     total_score += i
# Write the data to a text file
with open(output_file, 'w') as file:
    file.write('pept_res'+ '\t' + 'tar_res'+'\t'+'distance'+'\t'+'type'+'\t'+'SCORE'+'\n')
    for i in range(len(dict['distance'])):
        file.write(f"{dict['pept_res'][i]}\t{dict['tar_res'][i]}\t{dict['distance'][i]}\t{dict['type'][i]}\t{dict['score'][i]}\n")
    file.write(f'FINAL SCORE\t\t\t\t\t{total_score}')

#We append the scores
score_file='scores.csv'
# Check if the file exists
if not os.path.exists(score_file):
    # File doesn't exist, create it and write the header
    with open(score_file, 'w') as file:
        file.write('description'+ ',' + 'INTERACTING FILE'+ ',' + 'SCORE'+ '\n')

with open(score_file, 'a') as file:
     file.write(f"{protein_name},{args.csv},{total_score}\n")


#close pymol
pymol.cmd.quit()



