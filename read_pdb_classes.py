from Bio.PDB import *
import pandas as pd
import matplotlib.pyplot as plt 
from classes_definitions import *
import os
from Bio.PDB.SASA import ShrakeRupley
import sys
import pandas as pd
import fastparquet


class TrainingSet:
        def __init__(self, pdb_path):
            self.pdb_path = pdb_path

        def get_pdb_files(self):
            pdb_files = [(self.pdb_path + "/" + f) for f in os.listdir(self.pdb_path) if f.endswith('.pdb')]
            return pdb_files

        def get_formated_set(self):
            result_df = pd.DataFrame()
            for pdb in self.get_pdb_files():
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
            return result_df   
                        

