#!/usr/bin/env python

from meeko import MoleculePreparation
from rdkit import Chem
import sys

ligand_sdf_file = sys.argv[1]
out_file = sys.argv[2]

with Chem.SDMolSupplier(ligand_sdf_file) as r:
    molecule = next(r)
molecule = Chem.AddHs(molecule, addCoords=True)
meeko_prep = MoleculePreparation()
meeko_prep.prepare(molecule)
meeko_prep.write_pdbqt_file(out_file)
