#!/usr/bin/env python

"""
The functions uncharge(), run_rxn(), and protonator() are based on https://github.com/jensengroup/protonator .
"""

from pathlib import Path
from collections import Counter
import sys
import shutil
import re
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from io import StringIO
import csv

naming = sys.argv[1]
ref_sdf_files = sys.argv[2:]


def uncharge(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()

    return mol


def run_rxn(reactant, smarts):
    ps = True
    while ps:
        rxn = AllChem.ReactionFromSmarts(smarts)
        ps = rxn.RunReactants((reactant,))
        if ps:
            reactant = ps[0][0]
            Chem.SanitizeMol(reactant)

    return reactant


def protonator(m):
    rxns = [
        '[NX3;H2;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$(NCC[ND2]);!$(NC=N):1]>>[NH3+:1]',
        # Primary aliphatic amines
        '[NX3;H1;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$([N;R][C;R][O,S,P;R]);!$([N;R][C;R][C;R][O,S,P;R]);!$(NC[N+]);!$(NCC[N+]);!$(NC=N):1]>>[NH2+:1]',
        # Secondary aliphatic amines
        '[NX3;H0;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$([N;R][C;R][O,S,P;R]);!$([N;R][C;R][C;R][O,S,P;R]);!$(NC[N+]);!$(NCC[N+]);!$(NC=N):1]>>[NH+:1]',
        # Tertiary aliphatic amines
        '[NX2;H0;!$(NC=[O,S]);!$(Na);!$(*=[N,S,O,P]);!$(*=[C;R][S;R]):1]>>[NH+:1]',  # sp2 N
        '[NX2;H1;!$(NC=[O,S]);!$(Na);!$(*=[N,S,O,P]);!$(*=[C;R][S;R]):1]>>[NH2+:1]',  # sp2 N
        # '[nX3;H1;$(*(n)n):1]>>[n-;H0:1]', #1,2,3-triazole
        '[nX3;H1;$(*(n)nn):1]>>[n-;H0:1]',  # tetrazole
        '[O;H1;$(O[C,S]=O):1]>>[O-;H0:1]',  # carboxylic and sulfinic acid
        '[O;H1;$(O[C,c;R]=[C,c;R][C,c]=O):1]>>[O-;H0:1]', '[O;H1;$(Occc=O):1]>>[O-;H0:1]',  # Phenol near C=O
        '[O;H1;$(Oaaa[nX3+]):1]>>[O-;H0:1]', '[O;H1;$(Oaa[nX3+]):1]>>[O-;H0:1]', '[O;H1;$(Oa[nX3+]):1]>>[O-;H0:1]',
        '[O;H1;$(O[nX3+]):1]>>[O-;H0:1]',  # Phenol near pyridinium
        '[O;H1;$(Ocnnn):1]>>[O-;H0:1]', '[O;H1;$(Onnn):1]>>[O-;H0:1]',  # Phenol on 1,2,3-triazole
        '[O;H1;$(Oc(c[F,Cl,Br,I])c[F,Cl,Br,I]):1]>>[O-;H0:1]',  # Phenol flanked by halides
        '[O;H1;$(OP(=[O,S])):1]>>[O-;H0:1]']  # Phosphate O

    for rxn in rxns:
        molb = run_rxn(m, rxn)
        if molb:
            m = molb

    return m


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


warnings = []
no_uncharge = []
no_protonation = []
no_embedding = []
preparation_failed = []

for ref in ref_sdf_files:
    print("now processing: " + ref)

    rdBase.LogToPythonStderr()
    stderr = sys.stderr
    sio = sys.stderr = StringIO()

    if naming == "default":
        ligand = ref.split("__")[1].split('.sdf')[0].split("_")[0]
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
                mol = Chem.MolFromSmiles(smiles)
                try:
                    mol = uncharge(mol)
                except Exception as e:
                    print("Error: cannot uncharge molecule")
                    print("error: " + str(e))
                    no_uncharge.append(ligand)

                try:
                    mol = protonator(mol)
                    print("molecule protonated")
                except Exception as e:
                    print("'Error: cannot protonate molecule'")
                    print("error: " + str(e))
                    no_protonation.append(ligand)

                try:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                except Exception as e:
                    print("Error: cannot embed molecule")
                    print("error: " + str(e))
                    no_embedding.append(ligand)

                with Chem.SDWriter(ligand_sdf_name) as w:
                    w.write(mol)

                if ligand_sdf_resnum_name != ligand_sdf_name:
                    shutil.copy(sdf_file, resnum_sdf_file)
                print("ligand preparation done")
            else:
                print(ref + " ligand preparation failed")
                preparation_failed.append(ligand)
        else:
            print("sdf file already exists")
            if ligand_sdf_resnum_name != ligand_sdf_name:
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
                mol = Chem.MolFromSmiles(smiles)
                try:
                    mol = uncharge(mol)
                except Exception as e:
                    print("Error: cannot uncharge molecule")
                    print("error: " + str(e))
                    no_uncharge.append(ligand)

                try:
                    mol = protonator(mol)
                    print("molecule protonated")
                except Exception as e:
                    print("Error: cannot protonate molecule")
                    print("error: " + str(e))
                    no_protonation.append(ligand)

                try:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                except Exception as e:
                    print("Error: cannot embed molecule")
                    print("error: " + str(e))
                    no_embedding.append(ligand)

                with Chem.SDWriter(ligand_preped) as w:
                    w.write(mol)

                print("ligand preparation done")
            else:
                print(ref + " ligand preparation failed")
                preparation_failed.append(ligand)
        else:
            print("sdf file already exists")

    if len(sio.getvalue()) != 0:
        warnings.append([ligand, sio.getvalue().strip()])
    print("warnings/errors:")
    print(sio.getvalue())
    sys.stderr = stderr

# print ligand preparation overview to log file and save to problems file
set_prep_failed = set(preparation_failed)
set_no_uncharge = set(no_uncharge)
set_no_protonation = set(no_protonation)
set_no_embedding = set(no_embedding)

preparation_failed = list(set_prep_failed)
no_uncharge = list(set_no_uncharge)
no_protonation = list(set_no_protonation)
no_embedding = list(set_no_embedding)
warnings.sort()

preparation_failed.sort()
no_uncharge.sort()
no_protonation.sort()
no_embedding.sort()
warnings.sort()

with open('problems_ligand_prep.txt', 'w+') as f:
    f.write('Ligand preparation failed for: ')
    f.writelines(','.join(preparation_failed))
    f.write('\nUncharging failed for: ')
    f.writelines(','.join(no_uncharge))
    f.write('\nProtonation failed for: ')
    f.writelines(','.join(no_protonation))
    f.write('\nEmbedding failed for: ')
    f.writelines(','.join(no_embedding))
    f.write('\nWarnings: \n')
    wr = csv.writer(f)
    wr.writerows(warnings)

print('\n')
print('*********************************')
print('*  Ligand preparation overview  *')
print('*********************************')
print('\n')
print('Ligand preparation failed:')
print(*preparation_failed, sep='\n')
print('\n')
print('Uncharging failed:')
print(*no_uncharge, sep='\n')
print('\n')
print('Protonation failed:')
print(*no_protonation, sep='\n')
print('\n')
print('Embedding failed:')
print(*no_embedding, sep='\n')
print('\n')
print('Warnings:')
print(*warnings, sep='\n')