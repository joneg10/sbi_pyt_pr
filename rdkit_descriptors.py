from rdkit import Chem
from rdkit.Chem import Descriptors


# Create a molecule object from the string
mol = Chem.MolFromSmiles("C")

all = Descriptors.CalcMolDescriptors(mol)

charge = Descriptors.MaxAbsPartialCharge(mol) # Maximum absolute partial charge, maximum partial charge, minimum absolute partial charge and minimum partial charge.

densityM = Descriptors.FpDensityMorgan1(mol)  # Morgan1, 2 and 3 fingerprint density. (No se qué es esto).

estateindex = Descriptors.MaxEStateIndex(mol) # Maximum E-state index, minimum E-state index, range of E-state indices and sum of E-state indices.