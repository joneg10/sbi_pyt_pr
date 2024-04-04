#!/usr/bin/env python3

def main(): 
        
        import argparse
        from wisefrog.environment_classes import StructureAnalysis
        import pandas as pd
        import torch
        import torch.nn as nn
        import subprocess
        import tempfile
        import sys
        import os 
        import requests
        import time
        import pkg_resources



        # create class to raise exception when file is not PDB

        class NotPDBFile(Exception):

                def __init__(self, file):
                        self.file = file

                def __str__(self):
                        return f"{self.file} is not a PDB file and could not be parsed."


        # mean from the training test

        mean = torch.tensor([3.4002e+01, 2.5525e+00, 1.8183e-01, 2.1489e-01, 1.6584e-01, 1.5862e-01,
        1.5208e-02, 9.9080e-04, 6.9211e-02, 6.2484e-01, 4.7259e-03, 1.7122e-01,
        6.0211e-03, 1.7709e-01, 2.6683e-02, 6.7682e-03, 2.1887e-02, 1.1795e-01,
        1.3157e-01, 1.3263e-01, 1.3149e-01, 2.5896e-02, 1.1662e-02, 2.9237e-02,
        9.6419e-03, 1.2512e-01, 7.2629e-03, 2.4276e-03, 4.9125e-03, 4.5131e-02,
        1.5346e-02, 4.4691e-03, 4.4286e-03, 1.3372e-02, 1.5627e-02, 1.1347e-02,
        7.1361e-03, 5.4788e-03, 6.4970e-03, 7.1330e-03, 5.4508e-03, 5.7619e-03,
        3.6518e-03, 2.4850e-03, 1.9137e-03, 1.9424e-03, 1.8992e-03, 1.8830e-03,
        1.9270e-03, 1.3970e-04])

        # standard deviation from the training test 

        std = torch.tensor([2.1407e+01, 2.8410e+00, 9.9229e-02, 8.5333e-02, 4.4860e-02, 4.7455e-02,
        1.9999e-02, 2.3514e-02, 3.0145e-02, 8.9659e-02, 1.1212e-02, 4.5903e-02,
        1.3150e-02, 5.0138e-02, 2.6208e-02, 1.3145e-02, 2.3256e-02, 2.9567e-02,
        3.4257e-02, 3.7944e-02, 3.6928e-02, 2.4788e-02, 1.6897e-02, 2.7302e-02,
        1.6252e-02, 3.7580e-02, 1.4664e-02, 7.8930e-03, 1.2039e-02, 1.4760e-02,
        1.9918e-02, 1.2549e-02, 1.0726e-02, 1.8346e-02, 2.0205e-02, 1.7484e-02,
        1.4208e-02, 1.3194e-02, 1.3462e-02, 1.3644e-02, 1.3376e-02, 1.3112e-02,
        1.0064e-02, 9.1373e-03, 6.8294e-03, 7.0227e-03, 6.9380e-03, 6.6800e-03,
        6.8912e-03, 2.1643e-03])

        # create empty model to load it with the parameters from the trained model

        model = nn.Sequential(
        nn.Linear(50, 74),
        nn.ReLU(),
        nn.Linear(74, 50),
        nn.ReLU(),
        nn.Linear(50, 25),
        nn.ReLU(),
        nn.Linear(25, 1),
        nn.Sigmoid()
        )


        # load model


        neural_network_path = pkg_resources.resource_filename(__name__,"neural_network_2603_1988_pdbs_6.2A.pytorch")
        model.load_state_dict(torch.load(neural_network_path))


        # parse arguments
        parser = argparse.ArgumentParser(description='Process the input PDB that is provided by the user, and show the ligand binding residues, atoms, and/or 3D visualization with Chimera.')

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

                try: 
                        input_structure = StructureAnalysis(pdb_file = input_pdb)


                except IndexError:
                        raise NotPDBFile(input_pdb)
                
                except ValueError:
                        raise NotPDBFile(input_pdb)
                
                else:

                        environment_df = pd.DataFrame.from_dict(input_structure.get_input_environments(), orient='index')

                        # drop columns that are not in the training set

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
                                sys.stdout.write("You can use this selection in chimera to visualize the binding atoms:\n\n")
                                sys.stdout.write(atoms_to_select + "\n\n")

                                with tempfile.NamedTemporaryFile(mode='w', delete=True, suffix='.cmd', dir=".") as f:
                                        f.write(f'open {input_pdb}\n')  # Open your PDB file
                                        f.write(f'select {atoms_to_select}\n')  # Select the atom
                                        f.write('color red sel\n')  # Color the selected atom red
                                        f.write('surf\n')
                                        f.write('surf sel\n')  # Color the selected atom red
                                        f.write('surf')
                                        
                                        f.flush()

                                        # Run Chimera with the temporary file
                                        try:
                                                subprocess.Popen(['chimera', f.name])
                                                time.sleep(5)
                                                
                                        except FileNotFoundError:
                                                print("Chimera is not installed. Please install Chimera to visualize the binding atoms.")

        
        # if args.pdb_file is file, predict file 

        if args.pdb_file[0:4] == "web/":

                URL = "https://files.rcsb.org/download/" + args.pdb_file[4:8] + ".pdb"
                
                response = requests.get(URL)
             
                if response.status_code == 404:
                        print(f"Error: The PDB file {args.pdb_file[4:8]} was not found.")
                else:
                        sys.stdout.write(f"Downloading from {URL}\n")
                        
                        with open(args.pdb_file[4:8]+".pdb", "w") as structure:
                                structure.write(response.text)
                        
                        predict_binding(structure.name)
                        
                        os.remove(structure.name)
        




        elif os.path.isfile(args.pdb_file):
                predict_binding(args.pdb_file)

        else:
                for pdb_file in os.listdir(args.pdb_file):
                        if pdb_file[-4:] == ".pdb":
                                predict_binding(args.pdb_file + "/" + pdb_file)

