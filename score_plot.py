
import os
import pandas as pd

# Set the directory path where your folders are located
directory_path = '/emdata/cchacon/TFEB/PMPNN/TFEB_crop'

# Initialize an empty list to store DataFrames
dfs = []

# Iterate through each folder in the directory
for folder in os.listdir(directory_path):
    folder_path = os.path.join(directory_path, folder)

    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Iterate through each file in the folder
        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)

            # Check if the file is a .sc file
            if file.endswith('.sc'):
                # Read the .sc file into a DataFrame and append it to the list
                df = pd.read_csv(file_path, sep='\s+')
                dfs.append(df)

# Concatenate all DataFrames into one
final_dataframe = pd.concat(dfs, ignore_index=True)

# Display or use the final DataFrame as needed
print(final_dataframe)


