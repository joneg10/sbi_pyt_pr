from Bio.PDB import *


parser = PDBParser(PERMISSIVE = True, QUIET = True)
data = parser.get_structure('10gs', 'pdb_ids/10GS.pdb/pdb10gs.ent')

# Get all atoms in the structure
atoms = list(data.get_atoms())

# Create a NeighborSearch object
ns = NeighborSearch(atoms)

# List of standard amino acids
standard_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


# Initialize an empty list to store the ligand names


ligand_residues = []

for residue in data.get_residues():
    resname = residue.get_resname()
    if resname not in standard_residues and resname != 'HOH':
        ligand_residues.append(resname)

print(ligand_residues)


#############
## ligand binding site
#############
ligand_binding_site_atoms = []
# Iterate over the residues in the structure
for residue in data.get_residues():
    residue_name = residue.get_resname()
    # Check if the residue is a ligand
    if residue_name in ligand_names:
        # Iterate over the atoms in the ligand
        for atom in residue.get_atoms():
            # find all atoms within 8A of the ligand atom
            nearby = ns.search(atom.coord, 8)
            # filter atoms that belong to the protein
            nearby_protein_atoms = [a for a in nearby if a.get_parent().get_resname() in standard_residues]
            ligand_binding_site_atoms.extend([(nearby_atom.get_name(), nearby_atom.coord, nearby_atom.get_parent().get_resname()) for nearby_atom in nearby_protein_atoms])



#################
#### Environments
#################
            
            
environments = {}

for residue in data.get_residues():
    if residue.get_resname() not in ligand_residues and residue.get_resname() != 'HOH':  # Exclude water residues
        for selected_atom in residue.get_atoms():
            nearby_atoms = ns.search(selected_atom.coord, 8)  # Change the radius if needed

            # Filter nearby atoms that are not water and are not in residues of ligand_names
            nearby_atoms = [atom for atom in nearby_atoms if atom.get_parent().get_resname() != 'HOH' and atom.get_parent().get_resname() not in ligand_residues]
            
            # Initialize the atom counts and counter outside of the loop
            atom_counts = {}
            atom_counter = 0
            
            for nearby_atom in nearby_atoms:
                atom_counter += 1

                if nearby_atom.get_id() in atom_counts:
                    atom_counts[nearby_atom.get_id()] += 1
                else:
                    atom_counts[nearby_atom.get_id()] = 1

                if str("total" + nearby_atom.element) in atom_counts:
                    atom_counts["total" + nearby_atom.element] += 1
                else:
                    atom_counts["total" + nearby_atom.element] = 1

            # Add the atom counts to the environments dictionary
            atom_proportions = {atom_id: count / atom_counter for atom_id, count in atom_counts.items()}
            environments[(selected_atom.get_serial_number(), selected_atom.get_parent().get_resname())] = atom_proportions




















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



