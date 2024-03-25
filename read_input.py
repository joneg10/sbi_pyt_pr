import argparse
from copy_of_classes import *
import pandas as pd
import torch
import torch.nn as nn
import subprocess
import tempfile
import os
import sys

# updated 2024/03/25

mean = torch.tensor([3.4326e+01, 2.6050e+00, 1.8401e-01, 2.1455e-01, 1.6547e-01, 1.5861e-01,
        1.4536e-02, 6.8653e-04, 6.8937e-02, 6.2669e-01, 4.5766e-03, 1.7025e-01,
        5.9161e-03, 1.7665e-01, 2.6430e-02, 6.6565e-03, 2.2254e-02, 1.1848e-01,
        1.3249e-01, 1.3356e-01, 1.3265e-01, 2.5828e-02, 1.1521e-02, 2.9603e-02,
        9.4883e-03, 1.2613e-01, 7.1888e-03, 2.5005e-03, 4.7794e-03, 4.5757e-02,
        1.4937e-02, 4.2974e-03, 4.2952e-03, 1.2944e-02, 1.6003e-02, 1.0906e-02,
        6.7250e-03, 5.1890e-03, 6.1702e-03, 6.9537e-03, 5.1348e-03, 5.5339e-03,
        3.3744e-03, 2.2789e-03, 1.8582e-03, 1.8711e-03, 1.8232e-03, 1.8539e-03,
        1.8738e-03, 1.3164e-04])



std = torch.tensor([2.1607e+01, 2.8752e+00, 9.8401e-02, 8.3578e-02, 4.3584e-02, 4.6237e-02,
        1.9245e-02, 2.2545e-02, 2.9891e-02, 8.8816e-02, 1.0981e-02, 4.4675e-02,
        1.3010e-02, 4.9064e-02, 2.5966e-02, 1.2899e-02, 2.3266e-02, 2.9163e-02,
        3.3618e-02, 3.7120e-02, 3.6138e-02, 2.4512e-02, 1.6595e-02, 2.7098e-02,
        1.6052e-02, 3.6808e-02, 1.4388e-02, 7.9099e-03, 1.1460e-02, 1.4909e-02,
        1.9387e-02, 1.2239e-02, 1.0425e-02, 1.7721e-02, 2.0242e-02, 1.6955e-02,
        1.3498e-02, 1.2621e-02, 1.2962e-02, 1.3342e-02, 1.2756e-02, 1.2632e-02,
        9.4203e-03, 8.3321e-03, 6.6400e-03, 6.8027e-03, 6.6902e-03, 6.5515e-03,
        6.7056e-03, 2.0907e-03])

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

input_structure = StructureAnalysis(pdb_file=args.pdb_file)

environment_df = pd.DataFrame.from_dict(input_structure.get_input_environments(), orient='index')

columns_to_drop = [col for col in environment_df.columns if col not in StructureAnalysis.column_names and col != "residue"]

environment_df.drop(columns_to_drop, axis=1, inplace=True)

residue_col = environment_df.pop("residue")

X_to_predict = torch.tensor(environment_df.values, dtype=torch.float32)
X_to_predict_normalized = (X_to_predict - mean) / std

### predict

model = nn.Sequential(
    nn.Linear(49, 74),
    nn.ReLU(),
    nn.Linear(74, 49),
    nn.ReLU(),
    nn.Linear(49, 25),
    nn.ReLU(),
    nn.Linear(25, 1),
    nn.Sigmoid()
)

model.load_state_dict(torch.load("neural_network_2403_more_layers.pytorch"))



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
                #f.write('surf sel\n')  # Color the selected atom red
                
                f.flush()

                # Run Chimera with the temporary file
                subprocess.run(['chimera', f.name])
