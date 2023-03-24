#!/usr/bin/env python

"""
This script is based on parts of DiffDocks baseline_gnina.py
https://github.com/gcorso/DiffDock/blob/main/baselines/baseline_gnina.py
as of 28.02.2023 with additional pdbqt preparation using meeko
"""

from meeko import MoleculePreparation
from rdkit import Chem
import sys

ligand_sdf_file = sys.argv[1]
out_file = sys.argv[2]

try:
    suppl = Chem.SDMolSupplier(sdf_file)
    molecule = suppl[0]
    Chem.SanitizeMol(mol)
    molecule = Chem.AddHs(molecule, addCoords=True)

    molecule.RemoveAllConformers()
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol, ps)
    if id == -1:
        print('rdkit pos could not be generated without using random pos. using random pos now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(molecule, ps)
        AllChem.MMFFOptimizeMolecule(molecule, confId=0)

    meeko_prep = MoleculePreparation()
    meeko_prep.prepare(molecule)
    meeko_prep.write_pdbqt_file(out_file)

except Exception as e:
    print(e)
