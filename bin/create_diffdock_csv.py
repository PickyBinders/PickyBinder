#!/usr/bin/env python

from ligand_utility import diffdock_csv
from pathlib import Path
import sys

naming = sys.argv[1]
receptor_hs = sys.argv[2]
ref_sdf_files = sys.argv[3:]

for ref in ref_sdf_files:
    if naming == 'default':
        pdb_chain, ligand_resnum = ref.split("__")
        sdf_file_name = Path() / f"{ligand_resnum}"
        if receptor_hs == 'no':
            pdb_file_name = Path() / f"{pdb_chain}_Hs.pdb"
        else:
            pdb_file_name = Path() / f"{pdb_chain}.pdb"

        diffdock_csv(sdf_file_name, pdb_file_name)

    else:
        receptor = ref.split("_")[0]
        ligand = ref.split(".")[0]
        sdf_file_name = Path() / f"{ligand}_preped.sdf"
        if receptor_hs == 'no':
            pdb_file_name = Path() / f"{receptor}_receptor_Hs.pdb"
        else:
            pdb_file_name = Path() / f"{receptor}_receptor.pdb"

        diffdock_csv(sdf_file_name, pdb_file_name)
