import argparse
from classes_definitions import *
import pandas as pd
import torch
import torch.nn as nn
import subprocess
import tempfile
import os
import sys

mean = torch.tensor([3.4787e+01, 3.0443e+00, 1.9361e-01, 2.1253e-01, 1.6394e-01, 1.6309e-01,
        1.3352e-02, 2.7850e-02, 6.3293e-01, 4.4168e-03, 1.6767e-01, 3.2025e-02,
        7.1537e-02, 6.0628e-03, 1.8101e-01, 2.7267e-02, 1.0240e-02, 1.1805e-01,
        1.3070e-01, 1.3029e-01, 1.2897e-01, 1.2917e-01, 7.1875e-03, 2.4488e-02,
        4.0618e-02, 1.5385e-02, 5.7021e-03, 5.3287e-03, 1.2753e-02, 7.6583e-03,
        5.2297e-03, 1.1875e-02, 1.1624e-02, 5.0469e-03, 1.7303e-02, 1.8148e-03,
        4.7217e-03, 2.9069e-03, 1.8127e-03, 1.7826e-03, 1.7981e-03, 1.8821e-03,
        1.8380e-03, 6.1939e-03, 4.3010e-03, 6.3468e-03, 6.5647e-03, 2.8375e-03])

std = torch.tensor([2.4775e+01, 2.7933e+00, 7.6304e-02, 6.3380e-02, 3.0415e-02, 3.2761e-02,
        1.3394e-02, 1.9076e-02, 6.5999e-02, 8.4308e-03, 3.1463e-02, 2.0603e-02,
        2.1160e-02, 9.2975e-03, 3.5243e-02, 1.8131e-02, 1.1839e-02, 1.9518e-02,
        2.2387e-02, 2.2143e-02, 2.4054e-02, 2.0577e-02, 9.7518e-03, 1.7913e-02,
        1.3887e-02, 1.3184e-02, 8.6615e-03, 8.5087e-03, 1.2413e-02, 1.0304e-02,
        8.5216e-03, 1.1931e-02, 1.2604e-02, 8.4529e-03, 1.5614e-02, 5.8573e-03,
        8.6552e-03, 6.2257e-03, 4.7201e-03, 4.7085e-03, 4.7130e-03, 4.7791e-03,
        4.7531e-03, 9.3402e-03, 7.2882e-03, 9.4651e-03, 9.6712e-03, 6.3077e-03])

parser = argparse.ArgumentParser(description='Process input PDB')

parser.add_argument('-f', '--file',
                    dest="pdb_file",
                    action="store",
                    required=True,
                    help="Input PDB file")

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

input_structure = StructureAnalysis("../pdb_ids/1B42.pdb")

environment_df = pd.DataFrame.from_dict(input_structure.get_input_environments(), orient='index')

columns_to_drop = [col for col in environment_df.columns if col not in StructureAnalysis.column_names and col != "residue"]
environment_df.drop(columns_to_drop, axis=1, inplace=True)

residue_col = environment_df.pop("residue")

X_to_predict = torch.tensor(environment_df.values, dtype=torch.float32)
X_to_predict_normalized = (X_to_predict - mean) / std

### predict

model = nn.Sequential(
    nn.Linear(48, 72),
    nn.ReLU(),
    nn.Linear(72, 48),
    nn.ReLU(),
    nn.Linear(48, 1),
    nn.Sigmoid()
)


model.load_state_dict(torch.load("neural_network.pytorch"))



# prediction

prediction = model(X_to_predict_normalized) 

prediction_codes = pd.DataFrame({'code': environment_df.index,
                                'prediction': prediction.detach().numpy().flatten(),
                                'residue': residue_col})

residues_output = set(prediction_codes[prediction_codes["prediction"] > 0.5]["residue"].values)

residues_output = sorted(residues_output, key=lambda x: int(x[3:]))

if args.output_atoms:
        if args.output_file:
                with open(args.output_file, "w") as f:
                        f.write("ATOM".ljust(0) + "RESIDUE".rjust(16) + "\n\n")
                        for atom in prediction_codes[prediction_codes["prediction"] > 0.5]["code"].values:
                               f.write(atom.ljust(0) + prediction_codes[prediction_codes["code"] == atom]["residue"].values[0].rjust(20-len(atom))+ "\n")

        else:
                sys.stdout.write("ATOM".ljust(0) + "RESIDUE".rjust(16) + "\n")
                for atom in prediction_codes[prediction_codes["prediction"] > 0.5]["code"].values:
                        sys.stdout.write(atom.ljust(0) + prediction_codes[prediction_codes["code"] == atom]["residue"].values[0].rjust(20-len(atom))+ "\n")

else:
        if args.output_file:
                with open(args.output_file, "w") as f:
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
                f.write(f'open {args.pdb_file}\n')  # Open your PDB file
                f.write(f'select {atoms_to_select}\n')  # Select the atom
                f.write('color red sel\n')  # Color the selected atom red
                f.flush()

                # Run Chimera with the temporary file
                subprocess.run(['chimera', f.name])
