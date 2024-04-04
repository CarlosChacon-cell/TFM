import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='file to remove lines', required=True)
parser.add_argument('--mode', help='mode of the xvg ', required=True)
parser.add_argument('--folder', help='folder to save the images', required=True)
parser.add_argument('--i', help='folder to save the images', required=True)
args = parser.parse_args()

file = args.file
output_folder = args.folder  # Specify the folder to save the plots

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

if args.mode == 'G':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Gyrate', 'x', 'y', 'z'])
    df.plot(x='time(ps)', y='Gyrate')
    plt.savefig(os.path.join(output_folder, f'gyrate_plot_{args.i}.png'))  # Save plot as an image
    plt.close()

# RMSD
if args.mode == 'R':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ns)', 'RMSD'])
    df.plot(x='time(ps)', y='RMSD')
    plt.savefig(os.path.join(output_folder, 'RMSD_plot.png'))  # Save plot as an image
    plt.close()

# Potential
if args.mode == 'E':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Potential'])
    df.plot(x='time(ps)', y='Potential')
    plt.savefig(os.path.join(output_folder, 'potential_plot.png'))  # Save plot as an image
    plt.close()

# Pressure
if args.mode == 'P':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Pressure'])
    df.plot(x='time(ps)', y='Pressure')
    plt.savefig(os.path.join(output_folder, 'pressure_plot.png'))  # Save plot as an image
    plt.close()

# Density
if args.mode == 'D':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Density'])
    df.plot(x='time(ps)', y='Density')
    plt.savefig(os.path.join(output_folder, 'density_plot.png'))  # Save plot as an image
    plt.close()

if args.mode == 'T':
    df = pd.read_csv(file, sep='\s+', header=None, names=['time(ps)', 'Temperature'])
    df.plot(x='time(ps)', y='Temperature')
    plt.savefig(os.path.join(output_folder, 'temperature_plot.png'))  # Save plot as an image
    plt.close()

