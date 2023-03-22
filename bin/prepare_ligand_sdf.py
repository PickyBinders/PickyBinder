#!/usr/bin/env python

from ligand_utility import smiles_to_3d_mol
from pathlib import Path
import sys
import shutil
import re
from rdkit import Chem
from rdkit import rdBase
from io import StringIO

ref_sdf_files = sys.argv[1:]


def mol_and_smiles_from_file(sdf_file, mol2_file):
    """
    Convert a sdf or mol2 file to a RDKIT mol and further to a smiles string
    """
    mol = Chem.MolFromMolFile(sdf_file, sanitize=False)
    problem = False

    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol)
    except Exception as e:
        print("ligand preparation using sdf file failed")
        print("error: " + str(e))
        smiles = str(e)
        problem = True

    if problem and mol2_file.exists():
        print("try to use mol2 file for ligand preparation")
        mol = Chem.MolFromMol2File(ligand_mol2_name, sanitize=False)
        problem = False
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol)
            problem = False
        except Exception as e:
            print("ligand preparation using mol2 file failed")
            print("error: " + str(e))
            smiles = str(e)
            problem = True

    return mol, smiles, problem


for ref in ref_sdf_files:
    print("now processing: " + ref)

    rdBase.LogToPythonStderr()
    stderr = sys.stderr
    sio = sys.stderr = StringIO()

    pattern = re.compile("(__)")

    if bool(re.search(pattern, ref)):
        ligand = ref.split("__")[1].split("_")[0]
        ligand_sdf_name = ligand + ".sdf"
        ligand_sdf_resnum_name = ref.split("__")[1]
        ligand_full_name = ref.split(".sdf")[0]
        ligand_mol2_name = ligand_full_name + ".mol2"

        sdf_file = Path() / f"{ligand_sdf_name}"
        resnum_sdf_file = Path() / f"{ligand_sdf_resnum_name}"
        mol2_file = Path() / f"{ligand_mol2_name}"

        if not sdf_file.exists():
            mol, smiles, problem = mol_and_smiles_from_file(ref, mol2_file)

            if not problem:
                smiles_to_3d_mol(smiles, ligand_sdf_name)
                shutil.copy(sdf_file, resnum_sdf_file)
                print("ligand preparation done")
            else:
                print(ref + " ligand preparation failed")
        else:
            print("sdf file already exists")
            shutil.copy(sdf_file, resnum_sdf_file)

    else:
        ligand = ref.split(".sdf")[0]
        ligand_sdf_name = ref
        ligand_preped = ligand + "_preped.sdf"
        ligand_mol2_name = ligand + ".mol2"

        sdf_file = Path() / f"{ligand_preped}"
        mol2_file = Path() / f"{ligand_mol2_name}"

        if not sdf_file.exists():
            mol, smiles, problem = mol_and_smiles_from_file(ref, mol2_file)

            if not problem:
                smiles_to_3d_mol(smiles, ligand_preped)
                print("ligand preparation done")
            else:
                print(ref + " ligand preparation failed")
        else:
            print("sdf file already exists")

    print("warnings/errors:")
    print(sio.getvalue())
    sys.stderr = stderr
