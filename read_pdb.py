from Bio.PDB import *
from read_pdb_functions import *
import pandas as pd
import matplotlib.pyplot as plt 
from classes_definitions import *
import os
from Bio.PDB.SASA import ShrakeRupley
import sys
import pandas as pd
import fastparquet
import os

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

# take all pdbs from a folder
pdb_files = []
folder_path = 'pdb_ids/'

# 
for file in os.listdir(folder_path):
    if file.endswith('.pdb'):
        pdb_files.append(os.path.join(folder_path, file))

# get the environments and save them in a file
# file = 1
# Create an empty DataFrame to store the results
# result_df = pd.DataFrame()

# for pdb in pdb_files:
#     analysis = StructureAnalysis(pdb)
#     env = analysis.get_structure_environments()
#     df = pd.DataFrame.from_dict(env, orient='index')
    
#     # Add missing columns to df before reordering
#     missing_columns = [col for col in result_df.columns if col not in df.columns]
#     for column in missing_columns:
#         df[column] = np.nan
    
#     # Reorder the columns of the DataFrame to match the order of the first analysis
#     if result_df.empty:
#         result_df = df
#     else:
#         df = df[result_df.columns]
#         # Append the DataFrame to the result DataFrame
#         result_df = result_df._append(df)
    
#     if file == 1:
#         result_df.to_csv('../output.csv', mode='w', header=True, index=True)
#     else:    
#         result_df.to_csv('../output.csv', mode='a', header=False, index=True)
    
#     sys.stdout.write(f'File {pdb} processed\n ')
#     file += 1
    
result_df = pd.DataFrame()

for pdb in pdb_files:
    analysis = StructureAnalysis(pdb)
    env = analysis.get_structure_environments()
    df = pd.DataFrame.from_dict(env, orient='index')
    
    # Add missing columns to df before reordering
    missing_columns = [col for col in result_df.columns if col not in df.columns]
    for column in missing_columns:
        df[column] = np.nan
    
    # Reorder the columns of the DataFrame to match the order of the first analysis
    if result_df.empty:
        result_df = df
    else:
        df = df[result_df.columns]
        # Append the DataFrame to the result DataFrame
        result_df = result_df._append(df)
    
    output_file = '../output.parquet'
    if os.path.exists(output_file):
        # If the output file already exists, read it and append the new data
        existing_df = pd.read_parquet(output_file)
        result_df = existing_df._append(result_df)
    
    # Write the result DataFrame to the Parquet file
    result_df.to_parquet(output_file)
    
    sys.stdout.write(f'File {pdb} processed\n ')    






