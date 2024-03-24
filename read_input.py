import argparse
from copy_of_classes import *
import pandas as pd
import torch
import torch.nn as nn
import subprocess
import tempfile
import os
import sys

# updated 2024/03/24

mean = torch.tensor([3.4345e+01, 2.6304e+00, 1.8431e-01, 2.1410e-01, 1.6560e-01, 1.5844e-01,
        1.4571e-02, 7.1124e-04, 6.8492e-02, 6.2622e-01, 4.4856e-03, 1.7045e-01,
        5.8019e-03, 1.7658e-01, 2.6096e-02, 6.7702e-03, 2.2410e-02, 1.1857e-01,
        1.3285e-01, 1.3401e-01, 1.3307e-01, 2.5721e-02, 1.1329e-02, 2.9526e-02,
        9.3545e-03, 1.2638e-01, 7.1572e-03, 2.6325e-03, 4.7908e-03, 4.5739e-02,
        1.4794e-02, 4.1441e-03, 4.2293e-03, 1.2855e-02, 1.6048e-02, 1.0802e-02,
        6.7290e-03, 5.1975e-03, 6.0932e-03, 7.0248e-03, 5.1333e-03, 5.5515e-03,
        3.3838e-03, 2.1583e-03, 1.8184e-03, 1.8373e-03, 1.7867e-03, 1.8270e-03,
        1.8448e-03])



std = torch.tensor([2.2123e+01, 2.8879e+00, 9.8564e-02, 8.3486e-02, 4.3713e-02, 4.6114e-02,
        1.9261e-02, 2.2521e-02, 2.9928e-02, 8.8976e-02, 1.0895e-02, 4.4829e-02,
        1.2948e-02, 4.9173e-02, 2.5955e-02, 1.3113e-02, 2.3358e-02, 2.9197e-02,
        3.3629e-02, 3.7183e-02, 3.6178e-02, 2.4467e-02, 1.6550e-02, 2.7177e-02,
        1.6014e-02, 3.6727e-02, 1.4427e-02, 8.0930e-03, 1.1412e-02, 1.5712e-02,
        1.9400e-02, 1.2105e-02, 1.0375e-02, 1.7687e-02, 2.0215e-02, 1.6800e-02,
        1.3541e-02, 1.2690e-02, 1.2846e-02, 1.3404e-02, 1.2825e-02, 1.2726e-02,
        9.4475e-03, 8.0587e-03, 6.5647e-03, 6.7093e-03, 6.6030e-03, 6.5003e-03,
        6.6364e-03])

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
