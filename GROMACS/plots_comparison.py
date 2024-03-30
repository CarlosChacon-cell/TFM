import glob
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Setting up argparse
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Folder name pattern', required=False)
parser.add_argument('--mode', help='Mode' )
args = parser.parse_args()

# Extract folder name pattern from command line arguments
if args.folder:
    folder_name = args.folder
else: # Default pattern if not provided
    folder_name="PolG*"
# Create a list to store DataFrames
dataframes = []

# Get a list of folders matching the provided pattern
folder_list = glob.glob(folder_name)

# Iterate over each folder
if args.mode == 'R':
    for folder in folder_list:
        try:
            # Read RMSD data from each folder
            df = pd.read_csv(f'{folder}/rmsd.xvg', delim_whitespace=True, header=None, names=["Time", "RMSD"])
            # Plot RMSD data
            plt.plot(df["Time"], df["RMSD"], label=folder)
        except Exception as e:
            print(f"Error reading data from {folder}: {e}")
            continue
         
    # Add labels and legend
    plt.xlabel("Time (ps)")
    plt.ylabel("Backbone RMSD (A)")
    plt.title("RMSD")
    plt.legend()

    # Show the plot
    plt.show()
if args.mode== 'G':
    for folder in folder_list:
        try:
            # Read RMSD data from each folder
            df = pd.read_csv(f'{folder}/gyrate.xvg', delim_whitespace=True, header=None, names=["Time", "Gyration", "x", "y", "z"])
            # Plot RMSD data
            plt.plot(df["Time"], df["Gyration"], label=folder)
        except Exception as e:
            print(f"Error reading data from {folder}: {e}")
            continue
         
    # Add labels and legend
    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (A)")
    plt.title("Gyration")
    plt.legend()

    # Show the plot
    plt.show() 





#    dataframes.append(df)
