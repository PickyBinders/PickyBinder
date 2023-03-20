#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import sys
from Bio.PDB import PDBParser

tankbind_scripts = sys.argv[1]
sys.path.insert(0, tankbind_scripts)

from feature_utils import get_protein_feature, get_clean_res_list

complex = sys.argv[2]
p2rank_predictions = sys.argv[3]
box_size = sys.argv[4]
receptor_pdb = sys.argv[5]

with open(box_size, 'r') as infile:
    size = infile.read().strip()

    size_x = "size_x = " + size
    size_y = "size_y = " + size
    size_z = "size_z = " + size

df = pd.read_csv(p2rank_predictions)
df.columns = df.columns.str.strip()
df["name"] = df["name"].str.strip()

for row in range(len(df.index)):
    pocket = df.loc[row, "name"]
    pocket_file = Path() / f"{complex}_{pocket}.box"
    print(pocket_file)

    center_x = 'center_x = ' + str(df.loc[row, "center_x"])
    center_y = 'center_y = ' + str(df.loc[row, "center_y"])
    center_z = 'center_z = ' + str(df.loc[row, "center_z"])

    box = [center_x, center_y, center_z, size_x, size_y, size_z]

    with open(pocket_file, 'w+') as f:
        f.write('\n'.join(box))
        f.write("\n")

# add a box for the whole protein
"""
The definition of the box for the whole protein is based on tankbinds script. 
https://github.com/luwei0917/TankBind/blob/main/examples/construction_PDBbind_training_and_test_dataset.ipynb
as of 28.02.2023
"""

try:
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", receptor_pdb)
    res_list = get_clean_res_list(s.get_residues(), verbose=False, ensure_ca_exist=True)

    protein_dict = {}
    protein_dict[receptor_pdb] = get_protein_feature(res_list)

    protein_center = [str(a.round(3)) for a in protein_dict[receptor_pdb][0].mean(axis=0).numpy()]

    center_x = 'center_x = ' + str(protein_center[0])
    center_y = 'center_y = ' + str(protein_center[1])
    center_z = 'center_z = ' + str(protein_center[2])

    size = str(100)

    size_x = "size_x = " + size
    size_y = "size_y = " + size
    size_z = "size_z = " + size

    box = [center_x, center_y, center_z, size_x, size_y, size_z]

    pocket_file = Path() / f"{complex}_pocketProteinCenter.box"

    with open(pocket_file, 'w+') as f:
        f.write('\n'.join(box))
        f.write("\n")
except Exception as e:
    print("Docking box for protein center failed")
    print(e)
