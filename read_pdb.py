from Bio.PDB import *
from read_pdb_functions import *
import pandas as pd
import matplotlib.pyplot as plt 
from classes_definitions import *
import os
from Bio.PDB.SASA import ShrakeRupley


#### APPROACH 1: FUNCTIONS #####

parser = PDBParser(PERMISSIVE = True, QUIET = True)
data = parser.get_structure('10gs', '../pdb_ids/10gs.pdb')


# Get the ligand residues

ligand_residues = get_ligands_from_structure(data)


#############
## ligand binding site
#############

ligand_binding_site_atoms = get_ligand_binding_site_atoms(data)


#################
#### Environments
#################
            
            
get_structure_environments(data)

# create data frame
df = pd.DataFrame.from_dict(get_structure_environments(data), orient='index')

df.to_csv('output.csv', index=True)



# smoothing the values
smoothed_values = df['ligand_binding_site'].rolling(window=30).mean()

# binarize
df["ligand_binary"] = [1 if value > 0.7 else 0 for value in df["ligand_binding_site"]]
df["smoothed_values"] = smoothed_values
# plot
plt.plot(df.index, df['is_lbs'])
plt.xticks([])
plt.show(block=False)




#save list of atoms that are considered to be in the ligand binding site
filtered_indices = list(df[df["is_lbs"] == "Yes"].index)

with open("list_of_indices_yes_no.txt","w") as f:
    f.write(str(filtered_indices).replace('[','').replace(']','').replace('\'',''))






























# # Create a new PDB structure
# structure = Structure.Structure('Ligand Binding Site')

# # Create a new model
# model = Model.Model(0)

# # Add the model to the structure
# structure.add(model)

# # Create a new chain
# chain = Chain.Chain('A')

# # Add the chain to the model
# model.add(chain)

# # Iterate over the atoms in the ligand binding site and create a structure
# for i, (atom_name, atom_coord, residue_name) in enumerate(ligand_binding_site_atoms):
#     # Create a new residue
#     residue = Residue.Residue((' ', i+1, ' '), residue_name, ' ')
#     # Add the residue to the chain
#     chain.add(residue)
    
#     # Create a new Atom object
#     atom = Atom.Atom(atom_name, atom_coord, 0.0, 1.0, ' ', atom_name, i)
#     # Add the atom to the residue
#     residue.add(atom)

# # Create a PDBIO object
# pdb_io = PDBIO()
# # Set the structure to be written
# pdb_io.set_structure(structure)
# # Save the structure to a PDB file
# pdb_io.save('ligand_binding_site.pdb')






#### APPROACH 2: CLASSES #####

analysis = StructureAnalysis('../pdb_ids/10gs.pdb')
analysis.get_ligands_from_structure()
analysis.get_ligand_binding_site_atoms()
analysis.get_structure_environments()



# take all pdbs from a folder
pdb_files = []
folder_path = '../pdb_ids/'

# 
for file in os.listdir(folder_path):
    if file.endswith('.pdb'):
        pdb_files.append(os.path.join(folder_path, file))

# get the environments and save them in a file
for pdb in pdb_files:
    analysis = StructureAnalysis(pdb)
    env = analysis.get_structure_environments()
    df = pd.DataFrame.from_dict(env, orient='index')
    df.to_csv('output.csv', mode='a', header=False, index=True)


analysis = StructureAnalysis('../pdb_ids/10gs.pdb')

analysis.get_structure_environments()








