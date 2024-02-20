from Bio.PDB import *
import active_site


parser = PDBParser(PERMISSIVE = True, QUIET = True)
data = parser.get_structure('10gs', 'pdb_ids/10GS.pdb/pdb10gs.ent')

# Get all atoms in the structure
atoms = list(data.get_atoms())

# Create a NeighborSearch object
ns = NeighborSearch(atoms)

# Select an atom

for residue in data.get_residues():
    if residue.get_id()[0] == ' ' and residue.get_resname() != 'HOH':  # Exclude water residues
        for atom in residue.get_atoms():
            selected_atom = atom  # Change this to select a different atom

            # Find all atoms within 10A of the selected atom
            nearby_atoms = ns.search(selected_atom.coord, 10.0)  # Change the radius if needed

            # Check if any nearby atom is a ligand
            is_near_ligand = False
            for nearby_atom in nearby_atoms:
                nearby_residue = nearby_atom.get_parent()
                if nearby_residue.get_resname() != 'HOH' and nearby_residue.get_resname() != 'VWW':  # Exclude nearby water residues and specify the ligand name
                    is_near_ligand = True
                    break

            # Print "Yes" or "No" based on whether the atom is near a ligand
            if is_near_ligand:
                print("Yes")
            else:
                print("No")




# List of standard amino acids
standard_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Initialize an empty list to store the ligand names
ligand_names = []

# Iterate over the residues in the structure
for residue in data.get_residues():
    residue_name = residue.get_resname()
    # Check if the residue is not a standard amino acid or water
    if residue_name not in standard_residues and residue_name != 'HOH':
        ligand_names.append(residue_name)



# Iterate over the atoms in the structure
nearby_protein_atoms = []

for atom in atoms:
    residue = atom.get_parent()
    residue_name = residue.get_resname()
    # Check if the atom belongs to a ligand
    if residue_name in ligand_names:
        # Find all atoms within 10A of the ligand atom
        nearby = ns.search(atom.coord, 10.0)  # Change the radius if needed
        # Filter out atoms that belong to the protein
        nearby_protein = [a for a in nearby if a.get_parent().get_resname() in standard_residues]
        nearby_protein_atoms.extend(nearby_protein)

# Print the nearby protein atoms
for atom in nearby_protein_atoms:
    print(f"{atom.get_parent().get_resname()} {atom.get_name()} {atom.get_coord()}")

