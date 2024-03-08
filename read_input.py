from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
import atom_dict
import math
import argparse
import sys
from classes_definitions import *
import pandas as pd

parser = argparse.ArgumentParser(description='Process some pdb files.')

parser.add_argument('-i', '--input',
                    dest = "infile",
                    action = "store",
                    default = None,
                    help = "Input FASTA formatted file")


input_structure = StructureAnalysis(sys.argv[1])

df = pd.DataFrame.from_dict(input_structure.get_input_environments(), orient='index')






