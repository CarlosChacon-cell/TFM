from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import glob

parser=argparse.ArgumentParser()
parser.add_argument('--folder', '-f', help='Folder in which the fasta files are')
args = parser.parse_args()

def calculate_aa_composition(sequence,aa_counts, total_aa):
    for aa in sequence:
        if aa in aa_counts:
            aa_counts[aa] += 1
        else:
            aa_counts[aa] = 1
        total_aa += 1

    return aa_counts,total_aa

import matplotlib.pyplot as plt

def plot_aa_composition(aa_percentages, folder):
    labels = list(aa_percentages.keys())
    values = list(aa_percentages.values())
    title = f'Amino Acid Composition of {folder}'

    # Define custom colors for the bars
    colors = ['skyblue', 'lightgreen', 'lightcoral', 'orange', 'lightpink',
              'lightsalmon', 'lightblue', 'lightyellow', 'lightgrey', 'lavender']

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values, color=colors, edgecolor='black', linewidth=1.5)
    plt.xlabel('Amino Acid', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.title(title, fontsize=14)

    # Add gridlines for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Increase font size of tick labels
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # Optionally, add a legend
    # plt.legend(bars, labels, loc='upper right', fontsize=10)

    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    plt.show()


def process_fasta_files(fasta_files):
    aa_counts = {}
    total_aa=0
    for fasta_file in fasta_files:
        sequences = SeqIO.parse(fasta_file, 'fasta')
        for record in sequences:

            aa_counts,total_aa = calculate_aa_composition(record.seq,aa_counts,total_aa)
    aa_composition = {aa: (count / total_aa) * 100 for aa, count in aa_counts.items()}
    plot_aa_composition(aa_composition, args.folder)
    print(aa_counts)


if __name__ == "__main__":
    fasta_files = glob.glob(f"{args.folder}/*fasta")
    process_fasta_files(fasta_files)
