#!/usr/bin/env python

"""
This script is based on parts of DiffDocks baseline_gnina.py
https://github.com/gcorso/DiffDock/blob/main/baselines/baseline_gnina.py
as of 28.02.2023
"""

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, MolToPDBFile
from scipy.spatial.distance import cdist
import sys
import numpy as np

sdf_file = sys.argv[1]
ligand = sys.argv[2]
autobox_add = sys.argv[3]

try:
    suppl = Chem.SDMolSupplier(sdf_file)
    mol = suppl[0]
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    mol.RemoveAllConformers()
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol, ps)
    if id == -1:
        print('rdkit pos could not be generated without using random pos. using random pos now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol, confId=0)
    MolToPDBFile(mol, ligand + "_rdkit.pdb")

    rdkit_lig_pos = mol.GetConformer().GetPositions()

    diameter_pocket = np.max(cdist(rdkit_lig_pos, rdkit_lig_pos))
    size = np.round(diameter_pocket + int(autobox_add) * 2, 3)

    boxSize_file = Path() / f"{ligand}_boxSize.txt"

    with open(boxSize_file, 'w+') as f:
        f.write(str(size))
except Exception as e:
    print(ligand + ':')
    print('Box size could not be determined because of: ')
    print(e)
    print('The box size is set to 40Å for ' + ligand + '.')

    size = 40

    boxSize_file = Path() / f"{ligand}_boxSize.txt"

    with open(boxSize_file, 'w+') as f:
        f.write(str(size))
