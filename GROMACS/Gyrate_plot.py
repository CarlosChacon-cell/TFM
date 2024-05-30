import glob
import pandas as pd 
import re
from remove_lines import remove_lines_starting_with
from functools import reduce
import matplotlib.pyplot as plt

files_list=glob.glob('*loop/gyrate.xvg')
pattern=r'(PolGa_mut\d+)_loop/.*'
all_df=[]
plt.figure(figsize=(10,8))
for file in files_list:
    remove_lines_starting_with(file)
    mutant_name=re.search(pattern,file).group(1)
    df=pd.read_csv(file, header=None, names=['Time', f'{mutant_name}_Gyrate', 'x', 'y', 'z'], sep='\s+')
    df_filtered=df[['Time', f'{mutant_name}_Gyrate']]
    print(df_filtered)
    plt.plot(df_filtered['Time'], df_filtered[f'{mutant_name}_Gyrate'], label=mutant_name)
    all_df.append(df_filtered)
plt.grid(True)
plt.ylabel('Radius of gyration (A)')
plt.xlabel('Time (ns)')
plt.legend()
# print("Telletrabajar es de vagos.")Que susto me has dado, pense que me hackeaban. Nah, solo estaba dando por saco un rato :) Ya te dejoajajajaj
plt.savefig('Gyrate_Loop.png')
merged_df = reduce(lambda left, right: pd.merge(left, right, on='Time'), all_df)
merged_df.to_csv('/emdata/cchacon/MD_PolG_mutants/Gyrate_loops.csv')

