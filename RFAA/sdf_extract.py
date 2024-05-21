import argparse

def parse_molecule_data(file_path):
    namelist=[]

    with open(file_path, 'r') as file:
        content = file.read()

    # Split the content by the "$$$$" delimiter
    blocks = content.split('$$$$')

    molecules = []
    for block in blocks:
        # Strip any leading/trailing whitespace
        block = block.strip()
        if block:
            # Extract the molecule data up to the zinc_id
            lines = block.split('\n')
            zinc_id = ""
            for i, line in enumerate(lines):
                if line.startswith('>  <zinc_id>'):
                    zinc_id = lines[i + 1].strip()
                    namelist.append(lines[i + 1].strip()) 
                    break
            molecule_data = '\n'.join(lines[:i])
            if zinc_id:
                molecules.append(molecule_data + '\n')

    return molecules,namelist

# Path to your input file
parser=argparse.ArgumentParser()
parser.add_argument('--file', '-f')
args=parser.parse_args()
file_path = args.file

# Parse the file and get the list of molecules
molecules,comp_names = parse_molecule_data(file_path)
# Print each molecule data (or you can process further as needed)
for i in range(len(molecules)):
    with open(f'{comp_names[i]}.sdf', 'w') as file:
        file.write(f'\n{molecules[i]}')
    