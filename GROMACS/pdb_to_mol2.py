from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_to_mol2_with_hydrogens(pdb_filename, mol2_filename):
    # Read the PDB file
    mol = Chem.MolFromPDBFile(pdb_filename, removeHs=False)
    if mol is None:
        raise ValueError("Could not parse the PDB file.")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Compute 3D coordinates if needed (optional but recommended)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Write the MOL2 file
    with Chem.SDWriter(mol2_filename) as writer:
        writer.write(mol)

# Example usage
pdb_filename = 'ligand.pdb'
mol2_filename = 'output.mol2'
pdb_to_mol2_with_hydrogens(pdb_filename, mol2_filename)