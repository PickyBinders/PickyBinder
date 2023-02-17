#!/usr/bin/env python

from ligand_utility import mol_to_smiles, smiles_to_3d_mol
from pathlib import Path
import sys
import shutil
import re

ref_sdf_files = sys.argv[1:]

for ref in ref_sdf_files:
    pattern = re.compile("(__)")
    if bool(re.search(pattern, ref)):
        ligand = ref.split("__")[1].split("_")[0] + ".sdf"
        ligand_resnum = ref.split("__")[1]

        ligand_sdf = Path() / f"{ligand}"
        ligand_resnum_sdf = Path() / f"{ligand_resnum}"

        try:
            if not ligand_sdf.exists():
                smiles = mol_to_smiles(ref)
                smiles_to_3d_mol(smiles, ligand)
                shutil.copy(ligand_sdf, ligand_resnum_sdf)

            else:
                shutil.copy(ligand_sdf, ligand_resnum_sdf)

        except Exception as e:
            print(ligand, e)
    else:
        ligand = ref
        ligand_preped = ligand.split(".")[0] + "_preped.sdf"
        ligand_sdf = Path() / f"{ligand_preped}"
        try:
            if not ligand_sdf.exists():
                smiles = mol_to_smiles(ref)
                smiles_to_3d_mol(smiles, ligand_preped)
        except Exception as e:
            print(ligand, e)
