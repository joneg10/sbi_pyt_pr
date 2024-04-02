"""
create_training_set.py

This module contains the TrainingSet class, which is used to create a formatted training set from a directory of PDB files.

Classes:
    TrainingSet: Represents a training set created from a directory of PDB files.

Methods:
    get_pdb_files(): Returns a list of all PDB files in the directory.
    get_formated_set(): Do feature extraction and return a formatted training set in pandas DataFrame format.
"""

from Bio.PDB import *
import pandas as pd
import matplotlib.pyplot as plt 
from environment_classes import *
import os
from Bio.PDB.SASA import ShrakeRupley
import sys
import pandas as pd






class TrainingSet:
        """
    A class used to represent a training set created from a directory of PDB files.

    ...

    Attributes
    ----------
    pdb_path : str
    a string representing the path to the directory containing the PDB files

    Methods
    -------
    get_pdb_files():
    Returns a list of all PDB files in the directory.
    get_formated_set():
    Returns a DataFrame representing the formatted training set.
    """
            
        def __init__(self, pdb_path):
            """
            Constructs all the necessary attributes for the TrainingSet object.

            Parameters
            ----------
            pdb_path : str
                a string representing the path to the directory containing the PDB files
            """
            self.pdb_path = pdb_path




        def get_pdb_files(self):
            """
            Returns a list of all PDB files in the directory.

            Returns
            -------
            list
                a list of strings representing the paths to the PDB files
            """

            pdb_files = [(self.pdb_path + "/" + f) for f in os.listdir(self.pdb_path) if f.endswith('.pdb')]
            return pdb_files

        def get_formated_set(self):
            """
            Returns a DataFrame representing the formatted training set.

            The DataFrame is created by analyzing each PDB file in the directory using the StructureAnalysis class. 
            Each row in the DataFrame represents a different PDB file.

            Returns
            -------
            DataFrame
                a pandas DataFrame representing the formatted training set
            """
            result_df = pd.DataFrame()
            file_counter = 1


            for pdb in self.get_pdb_files():
                sys.stdout.write(f'Processing file {pdb}\n')
                analysis = StructureAnalysis(pdb)
                if len(analysis.get_ligands_from_structure()) > 0:
                    
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

                        sys.stdout.write(f'File {pdb} processed\n ')
                        sys.stdout.write(f"File {file_counter} of {len(self.get_pdb_files())}\n") 
                        file_counter += 1
                        
                    
            return result_df   
                        

