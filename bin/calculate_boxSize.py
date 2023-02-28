#!/usr/bin/env python

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

    MolToPDBFile(mol, ligand + "_rdkit.pdb")

    rdkit_lig_pos = mol.GetConformer().GetPositions()

    diameter_pocket = np.max(cdist(rdkit_lig_pos, rdkit_lig_pos))
    size = np.round(diameter_pocket + int(autobox_add) * 2, 3)

    boxSize_file = Path() / f"{ligand}_boxSize.txt"

    with open(boxSize_file, 'w+') as f:
        f.write(str(size))
except Exception as e:
    print(e)


