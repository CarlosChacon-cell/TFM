

import pandas as pd 

df=pd.read_excel('/home/cchacon/Downloads/Table_BM5.5.xlsx', header=0 , skiprows=[0,1,3,166,227])

pdb_names=[]
for value in df['Complex']:
    pdb_names.append(value[:4])
print(pdb_names)

filename='/emdata/cchacon/20240610_new_cutre/pdb_namaes.txt'

with open (filename, 'w') as file:
    for value in pdb_names:
        file.write(value+'\n')
