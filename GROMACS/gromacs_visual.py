'''
@chacon
This is a code used to plot the xvg files that you get from a GROMACS simulation.

Inputs:

--file: The xvg file to plot
--mode: The mode of the xvg file, whether it is a potential energy, temperature, etc...
--folder: Folder in which the plots must be saved (it can create it, the default is plots)

Outputs:

--A png file inside the selected folder

'''


import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np
from remove_lines import remove_lines_starting_with

#Function definition
def prettify_plot(ax):
    # Increase font sizes
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title(ax.get_title(), fontsize=14)
    ax.set_xlabel(ax.get_xlabel(), fontsize=12)
    ax.set_ylabel(ax.get_ylabel(), fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.6)




parser = argparse.ArgumentParser()
parser.add_argument('--file', help='file to plot', required=True)
parser.add_argument('--mode', help='mode of the xvg ', required=True)
parser.add_argument('--folder', help='folder to save the images', required=False, default='plots')
args = parser.parse_args()

remove_lines_starting_with(args.file)

file = args.file
output_folder = args.folder  # Specify the folder to save the plots

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

if args.mode == 'G':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Gyrate', 'x', 'y', 'z'])
    ax=df.plot(x='time(ps)', y='Gyrate', color='blue')
    prettify_plot(ax)
    plt.savefig(os.path.join(output_folder, f'gyrate_plot_{args.i}.png'))  # Save plot as an image
    plt.close()

# RMSD
if args.mode == 'R':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ns)', 'RMSD'])
    ax=df.plot(x='time(ps)', y='RMSD')
    prettify_plot(ax)
    plt.savefig(os.path.join(output_folder, 'RMSD_plot.png'))  # Save plot as an image
    plt.close()

# Potential
if args.mode == 'E':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Potential'])
    ax=df.plot(x='time(ps)', y='Potential')
    prettify_plot(ax)
    plt.savefig(os.path.join(output_folder, 'potential_plot.png'))  # Save plot as an image
    plt.close()

# Pressure
if args.mode == 'P':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Pressure'])
    ax=df.plot(x='time(ps)', y='Pressure')
    mean=df['Pressure'].mean()
    print('The mean pressure is: ', mean)
    ax.axhline(y=mean, color='r', linestyle='--', label=f'Mean Pressure: {mean:.2f}')
    prettify_plot(ax)
    plt.savefig(os.path.join(output_folder, 'pressure_plot.png'))  # Save plot as an image
    plt.close()

# Density
if args.mode == 'D':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Density'])
    ax=df.plot(x='time(ps)', y='Density')
    mean=df['Density'].mean()
    print('The mean Density is: ', mean)
    ax.axhline(y=mean, color='r', linestyle='--', label=f'Mean Density: {mean:.2f}')
    prettify_plot(ax)
    plt.savefig(os.path.join(output_folder, 'density_plot.png'))  # Save plot as an image
    plt.close()
#Temperature
if args.mode == 'T':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Temperature'])
    mean=df['Temperature'].mean()
    print('The mean temperature is: ', mean)
    ax=df.plot(x='time(ps)', y='Temperature')
    prettify_plot(ax)
    ax.axhline(y=mean, color='r', linestyle='--', label=f'Mean Temperature: {mean:.2f}')
    plt.legend()
    plt.savefig(os.path.join(output_folder, 'temperature_plot.png'))  # Save plot as an image
    plt.close()

