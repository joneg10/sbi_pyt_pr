#!/usr/bin/env python3

import argparse
from environment_classes import *
import pandas as pd
import torch
import torch.nn as nn
import subprocess
import tempfile
import sys
from algorithm import mean, std, model
import os 
import requests
import time
import pkg_resources

if __name__ == "__main__":

        parser = argparse.ArgumentParser(description='Process input PDB')

        parser.add_argument('-f', '--file',
                        dest="pdb_file",
                        action="store",
                        required=True,
                        help="Input files: can be an individual file, a directory containing multiple pdb files, or a PDB ID that will be downloaded from PDB using the format web/PDB_ID.")

        parser.add_argument('-c', '--chimera',
                        dest="open_chimera",
                        action="store_true",
                        help="Open Chimera with the output")

        parser.add_argument('-a', '--atoms',
                        dest="output_atoms",
                        action="store_true",
                        help="Show the prediction per atoms.")

        parser.add_argument('-o', '--output',
                        dest="output_file",
                        action="store",
                        help="Output file to write the results")

        args = parser.parse_args()


        # load model

        

        neural_network_path = pkg_resources.resource_filename(__name__,"neural_network_2603_1988_pdbs_6.2A.pytorch")
        model.load_state_dict(torch.load(neural_network_path))




        # function: predict lbs in input structure


        def predict_binding(input_pdb):
                """
                Predicts binding residues based on the input PDB file.

                Args:
                        input_pdb (str): Path to the input PDB file.

                Returns:
                        None
                """

                sys.stdout.write(f"========File: {input_pdb}========\n\n")

                input_structure = StructureAnalysis(pdb_file = input_pdb)

                environment_df = pd.DataFrame.from_dict(input_structure.get_input_environments(), orient='index')

                columns_to_drop = [col for col in environment_df.columns if col not in StructureAnalysis.column_names and col != "residue"]

                environment_df.drop(columns_to_drop, axis=1, inplace=True)

                residue_col = environment_df.pop("residue")

                X_to_predict = torch.tensor(environment_df.values, dtype=torch.float32)
                X_to_predict_normalized = (X_to_predict - mean) / std

                ### predict


                # prediction

                prediction = model(X_to_predict_normalized) 

                prediction_codes = pd.DataFrame({'code': environment_df.index,
                                                'prediction': prediction.detach().numpy().flatten(),
                                                'residue': residue_col})

                residues_output = set(prediction_codes[prediction_codes["prediction"] > 0.5]["residue"].values)

                residues_output = sorted(residues_output, key=lambda x: int(x[3:]))

                # output

                
                if args.output_atoms:
                        if args.output_file:
                                with open(args.output_file, "a") as f:
                                        f.write(f"========File: {input_pdb}========\n\n")
                                        f.write("ATOM".ljust(0) + "RESIDUE".rjust(16) + "\n\n")
                                        for atom in prediction_codes[prediction_codes["prediction"] > 0.4]["code"].values:
                                                f.write(atom.ljust(0) + prediction_codes[prediction_codes["code"] == atom]["residue"].values[0].rjust(20-len(atom))+ "\n")

                        else:
                                sys.stdout.write("ATOM".ljust(0) + "RESIDUE".rjust(16) + "\n")
                                for atom in prediction_codes[prediction_codes["prediction"] > 0.5]["code"].values:
                                        sys.stdout.write(atom.ljust(0) + prediction_codes[prediction_codes["code"] == atom]["residue"].values[0].rjust(20-len(atom))+ "\n")

                else:
                        if args.output_file:
                                with open(args.output_file, "a") as f:
                                        f.write(f"========File: {input_pdb}========\n\n")
                                        f.write(f"RESIDUE:\n\n")
                                        for residue in residues_output:
                                                f.write(residue+"\n")
                        else:
                                sys.stdout.write(f"RESIDUE:\n\n")
                                for residue in residues_output:
                                        sys.stdout.write(residue+"\n")


                atoms_to_select = str(list(prediction_codes[prediction_codes["prediction"] > 0.5]["code"].values)).replace("'", "").replace("[", "").replace("]", "").replace(" ", "")

                # Open Chimera
                if args.open_chimera:
                        with tempfile.NamedTemporaryFile(mode='w', delete=True, suffix='.cmd', dir=".") as f:
                                f.write(f'open {input_pdb}\n')  # Open your PDB file
                                f.write(f'select {atoms_to_select}\n')  # Select the atom
                                f.write('color red sel\n')  # Color the selected atom red
                                f.write('surf\n')
                                f.write('surf sel\n')  # Color the selected atom red
                                f.write('surf')
                                
                                f.flush()

                                # Run Chimera with the temporary file
                                subprocess.Popen(['chimera', f.name])
                                time.sleep(5)

        
        # if args.pdb_file is file, predict file 

        if args.pdb_file[0:4] == "web/":

                URL = "https://files.rcsb.org/download/" + args.pdb_file[4:8] + ".pdb"
                response = requests.get(URL)
                
                if response.status_code == 200:
                        
                        sys.stdout.write(f"Downloading from {URL}\n")
                        
                        with open(args.pdb_file[4:8]+".pdb", "w") as structure:
                                structure.write(response.text)
                        
                                predict_binding(structure.name)
                        
                else:
                        sys.stdout.write(f"The structure {args.pdb_file[4:8]} could not be retrieved from PDB\n")



        elif os.path.isfile(args.pdb_file):
                predict_binding(args.pdb_file)

        else:
                for pdb_file in os.listdir(args.pdb_file):
                        if pdb_file[-4:] == ".pdb":
                                predict_binding(args.pdb_file + "/" + pdb_file)

