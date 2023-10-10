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
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from io import StringIO
import csv


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
        '[NX3;H2;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$(NCC[ND2]);!$('
        'NC=N):1]>>[NH3+:1]',
        # Primary aliphatic amines
        '[NX3;H1;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$([N;R][C;R][O,S,'
        'P;R]);!$([N;R][C;R][C;R][O,S,P;R]);!$(NC[N+]);!$(NCC[N+]);!$(NC=N):1]>>[NH2+:1]',
        # Secondary aliphatic amines
        '[NX3;H0;!$(NC=[O,S]);!$(Na);!$(*[N,S,O,P]);!$(NC=[N;H0]);!$(NCC(F)F);!$(NCC#[C,N]);!$([N;R][C;R][O,S,'
        'P;R]);!$([N;R][C;R][C;R][O,S,P;R]);!$(NC[N+]);!$(NCC[N+]);!$(NC=N):1]>>[NH+:1]',
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


def protonate_and_optimize(smiles):
    problem = False
    uncharge_failed = False
    protonation_failed = False
    embedding_failed = False

    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        print("Error: smiles could not be converted to molecule object")
        print("error: " + str(e))
        problem = True

    if not problem:
        try:
            mol = uncharge(mol)
        except Exception as e:
            print("Error: cannot uncharge molecule")
            print("error: " + str(e))
            uncharge_failed = True

        try:
            mol = protonator(mol)
            print("molecule protonated")
        except Exception as e:
            print("Error: cannot protonate molecule")
            print("error: " + str(e))
            protonation_failed = True

        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol, confId=0)
        except Exception as e:
            print("Error: cannot embed molecule")
            print("error: " + str(e))
            embedding_failed = True

        return mol, uncharge_failed, protonation_failed, embedding_failed
    else:
        mol = False
        return mol, uncharge_failed, protonation_failed, embedding_failed


def write_mol_to_sdf(mol, outfile_name):
    try:
        with Chem.SDWriter(outfile_name) as w:
            w.write(mol)
    except Exception as e:
        print("Error: cannot write " + outfile_name)
        print("error: " + str(e))


def process_ligand_file(input_ligand):
    preparation_failed = False
    no_uncharge = False
    no_protonation = False
    no_embedding = False

    if input_ligand.endswith(".sdf"):
        mol = Chem.MolFromMolFile(input_ligand, sanitize=False)
        problem = False

        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol)
        except Exception as e:
            print("Ligand preparation using sdf file failed")
            print("error: " + str(e))
            smiles = str(e)
            problem = True

        if problem:
            mol2_file_name = input_ligand.replace(".sdf", ".mol2")
            mol2_file_path = Path() / f"{mol2_file_name}"

            if mol2_file_path.exists():
                print("Using mol2 file for ligand preparation")
                mol = Chem.MolFromMol2File(mol2_file_name, sanitize=False)
                try:
                    Chem.SanitizeMol(mol)
                    mol = Chem.RemoveHs(mol)
                    smiles = Chem.MolToSmiles(mol)
                    problem = False
                except Exception as e:
                    print("Ligand preparation using mol2 file failed")
                    print("error: " + str(e))
                    smiles = str(e)
                    problem = True

        if not problem:
            final_mol, no_uncharge, no_protonation, no_embedding = protonate_and_optimize(smiles)
            if final_mol is not False:
                write_mol_to_sdf(final_mol, ligand_preped)
            else:
                print("Ligand preparation failed for " + input_ligand)
                preparation_failed = True
        else:
            print("Ligand preparation failed for " + input_ligand)
            preparation_failed = True

    elif input_ligand.endswith(".mol2"):
        mol = Chem.MolFromMol2File(input_ligand, sanitize=False)
        problem = False

        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol)
        except Exception as e:
            print("Ligand preparation using mol2 file failed")
            print("error: " + str(e))
            smiles = str(e)
            problem = True

        if not problem:
            final_mol, no_uncharge, no_protonation, no_embedding = protonate_and_optimize(smiles)
            if final_mol is not False:
                write_mol_to_sdf(final_mol, ligand_preped)
            else:
                print("Ligand preparation failed for " + input_ligand)
                preparation_failed = True
        else:
            print("Ligand preparation failed for " + input_ligand)
            preparation_failed = True

    elif input_ligand.endswith(".smi"):
        problem = False
        try:
            with open(input_ligand) as f:
                line = f.readline()
                smiles = line.strip()
        except Exception as e:
            print("Failed to read smiles file")
            print("error: " + str(e))
            problem = True

        if not problem:
            final_mol, no_uncharge, no_protonation, no_embedding = protonate_and_optimize(smiles)
            if final_mol is not False:
                write_mol_to_sdf(final_mol, ligand_preped)
            else:
                print("Ligand preparation failed for " + input_ligand)
                preparation_failed = True

        else:
            print("Ligand preparation failed for " + input_ligand)
            preparation_failed = True

    else:
        smiles = input_ligand
        final_mol, no_uncharge, no_protonation, no_embedding = protonate_and_optimize(smiles)
        if final_mol is not False:
            write_mol_to_sdf(final_mol, ligand_preped)
        else:
            print("Ligand preparation failed for " + ligand_name)
            preparation_failed = True

    return preparation_failed, no_uncharge, no_protonation, no_embedding


if __name__ == "__main__":
    input_ligand = sys.argv[1]
    ligand_name = sys.argv[2]
    ligand_preped = ligand_name + "_preped.sdf"

    rdBase.LogToPythonStderr()
    stderr = sys.stderr
    sio = sys.stderr = StringIO()
    warnings = []

    try:
        preparation_failed, no_uncharge, no_protonation, no_embedding = process_ligand_file(input_ligand)
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

    # print warnings
    if len(sio.getvalue()) != 0:
        warnings.append([input_ligand, sio.getvalue().strip()])
    print("warnings/errors:")
    print(sio.getvalue())
    sys.stderr = stderr

    # print ligand preparation overview to log file
    print("\n")
    print("*********************************")
    print("*  Ligand preparation overview  *")
    print("*********************************")
    print("\n")
    print("Ligand preparation failed:")
    if preparation_failed:
        print(input_ligand)
    print("\n")
    print("Uncharging failed:")
    if no_uncharge:
        print(input_ligand)
    print("\n")
    print("Protonation failed:")
    if no_protonation:
        print(input_ligand)
    print("\n")
    print("Embedding failed:")
    if no_embedding:
        print(input_ligand)
    print("\n")
    print("Warnings:")
    print(*warnings, sep="\n")
