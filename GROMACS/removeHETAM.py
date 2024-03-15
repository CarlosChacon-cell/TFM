

#REMOVE DDS 
import argparse 


parser = argparse.ArgumentParser()
parser.add_argument('--file', help='file to remove lines', required=True)
parser.add_argument('--tag', help='mode of the xvg ', required=True)
args = parser.parse_args()

def remove_unwanted_lines(input_pdb_file):
    try:
        # Read the contents of the input PDB file
        with open(input_pdb_file, 'r') as f:
            pdb_lines = f.readlines()

        # Filter out lines with residue name DDS
        pdb_lines_filtered = [line for line in pdb_lines if not (line.startswith('HETATM') and line[17:20].strip() == args.tag)]

        # Write the filtered lines to a new PDB file
        output_pdb_file = input_pdb_file.split('.')[0] + "_modified.pdb"
        with open(output_pdb_file, 'w') as f:
            f.writelines(pdb_lines_filtered)

        print("Modified PDB file saved as:", output_pdb_file)
    except FileNotFoundError:
        print("Error: Input file not found.")
    except Exception as e:
        print("An error occurred:", e)

def main():
    input_pdb_file = args.file
    remove_unwanted_lines(input_pdb_file)

if __name__ == "__main__":
    main()
