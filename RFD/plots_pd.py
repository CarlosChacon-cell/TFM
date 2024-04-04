#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import shutil
import numpy as np
import argparse
import re
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--extract_hits", type=int, help = "Use 1 for extracting hits")
parser.add_argument("--folder", help = "Folder names, use an * and use the same structure for every folder name !!")
parser.add_argument("--pd", type=int, help = "output from partial diffusion use 1 or 0")
args, unknown = parser.parse_known_args()



# Columns to use for the scatter plot
y_column = 'plddt_binder'
x_column = 'pae_interaction'

# Initialize an empty list to store dataframes
dataframes = []

folder_list=glob.glob(args.folder)
print(folder_list)
#store input values

targets=['pae_interaction', 'plddt_binder']
extracted_values={'pae_interaction':[],'plddt_binder':[]}


# # Concatenate all dataframes into a single dataframe
for i in folder_list:
    file_list = glob.glob(f'{i}/output/run_*/*.sc', recursive=True)
    for file_path in file_list:
        with open(file_path, 'rb') as input_file:
            df = pd.read_table(file_path, sep=r'\s+')
            df = df[df != 'plddt_binder'].dropna()
            df['campaign'] = i
            df['len'] = 0
            for index, row in df.iterrows():
                run_number = re.search(r'run_\d+', row['description']).group()
                descriptor= row['description'].split('_dldesign')[0]
                file_name = f'{i}/output/'+f'{run_number}'+'/'+f'{descriptor}'+'.pdb'
                length = subprocess.check_output("grep \' CA \' "+file_name+" |  wc -l", shell=True)
                length = float(length) - float(575)
                df.at[index, 'len'] = length
            #print(df)
            dataframes.append(df)
    
    if args.pd==1:
        input_set=glob.glob(f'{i}/inputs/*.pdb')
        for inp in input_set:
            with open (inp, 'r') as protein:
                for line in protein:
                    for target in targets:
                        if line.startswith(target):
                            value = float(line.split()[1])
                            extracted_values[target].append(value)

        print(extracted_values)


combined_df = pd.concat(dataframes, ignore_index=True)
combined_df.to_csv('output_af2.csv', index=False)
print(combined_df)

#Define colors for different values in the 'Hotspot' column
color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:gray', 'tab:pink', 'tab:olive']
colors = dict(zip(folder_list, color_list[:len(folder_list)]))


combined_df['pae_interaction'] = combined_df['pae_interaction'].astype(float)
combined_df['plddt_binder'] = combined_df['plddt_binder'].astype(float)



if args.extract_hits: extract_hits = True
else: extract_hits = False

if extract_hits: 
    print('Hits will be extracted')
    hits = combined_df[(combined_df['pae_interaction'] <= 18) & (combined_df['plddt_binder'] >= 80)]
    for index, row in hits.iterrows():
        print('###############################\nEXTRACTING HIT\n###############################\n' + str(row))
        run_number = re.search(r'run_\d+', row['description']).group()
        os.system('silentextractspecific'+' '+row['campaign']+'/output/'+ run_number + '/' + run_number + '_input_out_af2.silent ' + row['description'] + ' > extraction.log')

elif extract_hits == 0:
    print('Hits until now (you can extract hits by using --extract_hits 1):')
    hits = combined_df[(combined_df['pae_interaction'] <= 10) & (combined_df['plddt_binder'] >= 80)]
    print('Filename\tCampaign')
    for index, row in hits.iterrows():
        print(str(row['description']) + '\t' + str(row['campaign']))


# Create a scatterplot
# plt.figure(figsize=(8, 6))
        
for campaign, color in colors.items():
    print(campaign, color)
    subset = combined_df[combined_df['campaign'] == campaign]
    plt.scatter(subset['pae_interaction'], subset['plddt_binder'], label=campaign, c=color)

# Add labels and legend
plt.xlabel('pae_interaction')
plt.ylabel('plddt_binder')
plt.title(f'Total number: {len(combined_df)}')
plt.legend(title='Campaign')
plt.gca().invert_xaxis()
plt.axhline(y=80, color='gray', linestyle='--')
plt.axvline(x=10, color='gray', linestyle='--')

# Show the plot
plt.show()

# Calculate the number of rows and columns based on the number of unique hotspots
num_unique_hotspots = len(colors)
num_cols = 3
num_rows = (num_unique_hotspots + num_cols - 1) // num_cols

# Create a single figure with subplots
fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows), gridspec_kw={'hspace': 0.5})


# Flatten the axs array to make it easier to iterate over
axs = axs.ravel()

#We create a counter to point out the inputs values
counter=0

# Iterate through unique values in 'Hotspot' and create scatterplots
for i, (hotspot, color) in enumerate(colors.items()):
    subset = combined_df[combined_df['campaign'] == hotspot]
    axs[i].scatter(subset['pae_interaction'], subset['plddt_binder'], c=color)
    axs[i].set_title(f'{hotspot} - Total: {len(subset)}')
    axs[i].set_xlabel('pae_interaction')
    axs[i].set_ylabel('plddt_binder')
    axs[i].invert_xaxis()
    axs[i].axhline(y=80, color='gray', linestyle='--')
    axs[i].axvline(x=10, color='gray', linestyle='--')
    if args.pd==1:
        axs[i].scatter(extracted_values['pae_interaction'][counter], extracted_values['plddt_binder'][counter], c='black')
        counter=counter+1

# Hide any empty subplots
for i in range(num_unique_hotspots, num_rows * num_cols):
    fig.delaxes(axs[i])

# Show the plot
plt.show()


# Calculate the number of rows and columns needed for subplots
num_rows = len(colors)
num_cols = 2  # You can change this if you want more columns


# Create a single figure with subplots
fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 6 * num_rows), gridspec_kw={'hspace': 0.5})

for i, (hotspot, color) in enumerate(colors.items()):
    subset = combined_df[combined_df['campaign'] == hotspot]
    axs[i,0].hist(subset[x_column], color=color, bins = 30, density = True )
    axs[i,0].set_title(f'{hotspot} - Total: {len(subset)}')
    axs[i,0].set_xlabel(x_column)
    axs[i,0].axvline(x=10, color='gray', linestyle='--')
    axs[i,0].set_ylabel('Density')

for i, (hotspot, color) in enumerate(colors.items()):
    subset = combined_df[combined_df['campaign'] == hotspot]
    axs[i,1].hist(subset[y_column], color=color, bins = 30, density = True )
    axs[i,1].set_xlabel(y_column)
    axs[i,1].set_title(f'{hotspot} - Total: {len(subset)}')
    axs[i,1].axvline(x=80, color='gray', linestyle='--')
    axs[i,1].set_ylabel('Density')

# Show the plot
plt.show()


if args.pd != 1:
    # Create a single figure with subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows), gridspec_kw={'hspace': 0.5})


    # Flatten the axs array to make it easier to iterate over
    axs = axs.ravel()


    # Iterate through unique values in 'Hotspot' and create scatterplots
    for i, (hotspot, color) in enumerate(colors.items()):
        subset = combined_df[combined_df['campaign'] == hotspot]
        axs[i].scatter(subset['pae_interaction'], subset['len'], c=color)
        axs[i].set_title(f'{hotspot} - Total: {len(subset)}')
        axs[i].set_xlabel('pae_interaction')
        axs[i].set_ylabel('length')
        axs[i].axvline(x=10, color='gray', linestyle='--')

    # Hide any empty subplots
    for i in range(num_unique_hotspots, num_rows * num_cols):
        fig.delaxes(axs[i])

    # Show the plot
    plt.show()

