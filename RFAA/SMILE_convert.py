from rdkit import Chem
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--file', '-f' , help='file to convert')
args=parser.parse_args()

mol_path=args.file
print(mol_path)
# Read PDBQT file
pdbqt_supplier = Chem.rdmolfiles.ForwardSDMolSupplier(mol_path, removeHs=False, sanitize=False)
mol = next(pdbqt_supplier)

mol_name=mol_path.split('.')[0]
# Convert to SMILES
smiles = Chem.MolToSmiles(mol)
with open(f'{mol_name}.smiles', 'w') as f:
    f.write(smiles)

# Convert to SDF
w = Chem.SDWriter(f'{mol_name}.sdf')
w.write(mol)
w.close()
