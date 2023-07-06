#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import sys
from Bio.PDB import PDBParser

tankbind_scripts = sys.argv[1]
sys.path.insert(0, tankbind_scripts)

from feature_utils import get_protein_feature, get_clean_res_list


def whole_protein_box(complex_name, pdb_file):
    """
    The definition for the box using the whole protein is based on tankbinds script.
    https://github.com/luwei0917/TankBind/blob/main/examples/construction_PDBbind_training_and_test_dataset.ipynb
    as of 28.02.2023

    Returns a text file containing the protein's center coordinates and the box diameter (100Å) that can
    be utilized for docking with Autodock Vina-like tools.
    """
    try:
        parser = PDBParser(QUIET=True)
        s = parser.get_structure("x", pdb_file)
        res_list = get_clean_res_list(s.get_residues(), verbose=False, ensure_ca_exist=True)

        protein_dict = {pdb_file: get_protein_feature(res_list)}

        protein_center = [str(a.round(3)) for a in protein_dict[pdb_file][0].mean(axis=0).numpy()]

        center_x = 'center_x = ' + str(protein_center[0])
        center_y = 'center_y = ' + str(protein_center[1])
        center_z = 'center_z = ' + str(protein_center[2])

        size = str(100)

        size_x = "size_x = " + size
        size_y = "size_y = " + size
        size_z = "size_z = " + size

        box = [center_x, center_y, center_z, size_x, size_y, size_z]

        pocket_file = Path() / f"{complex_name}_pocketProteinCenter.box"

        with open(pocket_file, 'w+') as f:
            f.write('\n'.join(box))
            f.write("\n")
    except Exception as e:
        print("Docking box for protein center failed")
        print(e)


def get_box_size(box_size_file):
    """
    Reads the box_size file and returns the size.
    """
    with open(box_size_file, 'r') as infile:
        size = infile.read().strip()

    return size


def box_per_predicted_pocket(complex_name, box_size, p2rank_predictions):
    """
    Writes for each P2Rank predicted pocket a docking box file with the given box size.

    Parameters:
        complex_name (str):        name of the complex, used for the out file name
        box_size (int):            size of the box diameter
        p2rank_predictions (file): prediction.csv from P2Rank

    Returns:
        pocket_file: contains the predicted pocket center coordinates and the box size in the form of newline separated
                     center_x = , center_y = , center_z = , size_x = , size_y = , size_z = , which can be utilized
                     for docking with Autodock Vina-like tools.
    """
    df = pd.read_csv(p2rank_predictions)
    df.columns = df.columns.str.strip()
    df["name"] = df["name"].str.strip()

    for row in range(len(df.index)):
        pocket = df.loc[row, "name"]
        pocket_file = Path() / f"{complex_name}_{pocket}.box"
        print(pocket_file)

        center_x = 'center_x = ' + str(df.loc[row, "center_x"])
        center_y = 'center_y = ' + str(df.loc[row, "center_y"])
        center_z = 'center_z = ' + str(df.loc[row, "center_z"])

        size_x = "size_x = " + box_size
        size_y = "size_y = " + box_size
        size_z = "size_z = " + box_size

        box = [center_x, center_y, center_z, size_x, size_y, size_z]

        with open(pocket_file, 'w+') as f:
            f.write('\n'.join(box))
            f.write("\n")


def box_for_defined_bs(complex_name, box_size, bs_coordinates):
    """
    Writes a docking box file with the given binding site coordinates and box size.

    Parameters:
        complex_name (str):   name of the complex, used for the out file name
        box_size (int):       size of the box diameter
        bs_coordinates (str): 'center-x_center-y_center-z' (underscore-separated coordinates of the binding site center)

    Returns:
        pocket_file: contains the provided binding site coordinates and box size in the form of newline separated
                     center_x = , center_y = , center_z = , size_x = , size_y = , size_z = , which can be utilized
                     for docking with Autodock Vina-like tools.
    """
    center_x = 'center_x = ' + bs_coordinates.split('_')[0]
    center_y = 'center_y = ' + bs_coordinates.split('_')[1]
    center_z = 'center_z = ' + bs_coordinates.split('_')[2]

    size_x = "size_x = " + box_size
    size_y = "size_y = " + box_size
    size_z = "size_z = " + box_size

    box = [center_x, center_y, center_z, size_x, size_y, size_z]

    pocket = "pocketBS"
    pocket_file = Path() / f"{complex_name}_{pocket}.box"

    with open(pocket_file, 'w+') as f:
        f.write('\n'.join(box))
        f.write("\n")


def main():
    """
    This script writes docking box files for either binding pockets predicted using P2Rank and the whole protein or
    a defined binding site and the whole protein.

    Example P2Rank predictions:
      define_docking_boxes.py path_to_tankbind_scripts complex_name pdb_file box_size_file p2rank_predictions_csv false

    Example binding site:
      define_docking_boxes.py path_to_tankbind_scripts complex_name pdb_file box_size_file false bs_coordinates
    """
    complex = sys.argv[2]
    receptor_pdb = sys.argv[3]
    box_size_file = sys.argv[4]
    p2rank_predictions = sys.argv[5]
    coordinates = sys.argv[6]

    whole_protein_box(complex, receptor_pdb)

    box_size = get_box_size(box_size_file)

    if coordinates == 'false':
        box_per_predicted_pocket(complex, box_size, p2rank_predictions)
    else:
        box_for_defined_bs(complex, box_size, coordinates)


if __name__ == "__main__":
    main()
