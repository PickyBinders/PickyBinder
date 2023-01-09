#!/usr/bin/env python

from ligand_utility import mol_to_smiles, smiles_to_3d_mol
from pathlib import Path
import sys

ref_sdf_files = sys.argv[1:]

for ref in ref_sdf_files:
  
  #ligand = ref.split("__")[1]
  ligand = ref.split("__")[1].split("_")[0] + ".sdf"
  
  ligand_sdf = Path() / f"{ligand}"
  
  if not ligand_sdf.exists():
        smiles = mol_to_smiles(ref)
        smiles_to_3d_mol(smiles, ligand)
