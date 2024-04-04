
import numpy as np
# Read the file
with open("pae.txt", "r") as file:
    data = file.readlines()

# Initialize lists to store the values
all_values = []
values=[]

# Iterate through the lines in the file
for line in data:
    # Check if the line is not empty
    if line:
        for element in line.strip().split():
            if element != ']':
                try:
                    element=float(element.strip('['))
                    values.append(element)
                except ValueError:
                    if len(values) != 0:
                        all_values.append(values)
                        values=[]
                        continue
                    else:
                        continue
            else:
                continue

residues=[]
with open('residues.txt', 'r') as resfile:
    for line in resfile:
        residues.append(int(line))

print(residues)

filtered_pae=[]

for residue in residues:
    filtered_pae.append(all_values[residue])

mean=[]
for lists in filtered_pae:
    mean.append(np.mean(lists))
mean=np.mean(mean)

print('\n#############\n')
print('\n#############\n')

print(' MEAN PAE: ', mean)

print('\n#############\n')
print('\n#############\n')
