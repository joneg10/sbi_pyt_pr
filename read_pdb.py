from Bio.PDB import *


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
            nearby_atoms = ns.search(selected_atom.coord, 2)  # Change the radius if needed

            # Check if any nearby atom is a ligand
            is_near_ligand = False
            for nearby_atom in nearby_atoms:
                nearby_residue = nearby_atom.get_parent()
                if nearby_residue.get_resname() == 'VWW':  # Exclude nearby water residues and specify the ligand name
                    is_near_ligand = True
                    break

            # Print "Yes" or "No" based on whether the atom is near a ligand
            if is_near_ligand:
                continue
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



ligand_names

# Iterate over the atoms in the structure
ligand_binding_site_atoms = []





ligand_binding_site_atoms = []
# Iterate over the residues in the structure
for residue in data.get_residues():
    residue_name = residue.get_resname()
    # Check if the residue is a ligand
    if residue_name in ligand_names:
        # Iterate over the atoms in the ligand
        for atom in residue.get_atoms():
            # find all atoms within 5A of the ligand atom
            nearby = ns.search(atom.coord, 5)
            # filter atoms that belong to the protein
            nearby_protein_atoms = [a for a in nearby if a.get_parent().get_resname() in standard_residues]
            ligand_binding_site_atoms.extend([(nearby_atom.get_name(), nearby_atom.coord, nearby_atom.get_parent().get_resname()) for nearby_atom in nearby_protein_atoms])




# Create a new PDB structure
structure = Structure.Structure('Ligand Binding Site')

# Create a new model
model = Model.Model(0)

# Add the model to the structure
structure.add(model)

# Create a new chain
chain = Chain.Chain('A')

# Add the chain to the model
model.add(chain)



# Iterate over the atoms in the ligand binding site and create a structure
for i, (atom_name, atom_coord, residue_name) in enumerate(ligand_binding_site_atoms):
    # Create a new residue
    residue = Residue.Residue((' ', i+1, ' '), residue_name, ' ')
    # Add the residue to the chain
    chain.add(residue)
    
    # Create a new Atom object
    atom = Atom.Atom(atom_name, atom_coord, 0.0, 1.0, ' ', atom_name, i)
    # Add the atom to the residue
    residue.add(atom)

# Create a PDBIO object
pdb_io = PDBIO()
# Set the structure to be written
pdb_io.set_structure(structure)
# Save the structure to a PDB file
pdb_io.save('ligand_binding_site.pdb')



