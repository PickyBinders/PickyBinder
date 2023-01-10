#!/usr/bin/env python

from pathlib import Path

import requests
import prody as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import RemoveStereochemistry, AssignStereochemistryFrom3D
from io import StringIO
import csv


def get_ligand_sdf(row, output_folder):
    pdb_chain, ligand_resnum = row.PocketID.split(":")
    ligand_resnum = ligand_resnum.replace(".", "_")
    ligand_link = f"https://models.rcsb.org/v1/{row.entry_pdbid.lower()}/ligand?auth_seq_id={row.resnum}&label_asym_id={row.said}&encoding=sdf&filename={row.entry_pdbid.lower()}_{row.said}_{row.Ligand}.sdf"
    request = requests.get(ligand_link, '\n')
    sdf_filename = Path(output_folder) / f"{pdb_chain}__{ligand_resnum}.sdf"
    with open(sdf_filename, 'w+') as f:
        f.write(request.text)


def get_protein_pdb(row, output_folder):
    pdb_chain, ligand_resnum = row.PocketID.split(":")
    pdb_file = Path(output_folder) / f"{row.entry_pdbid.lower()}.pdb"
    if not pdb_file.exists():
        protein_link = f"https://files.rcsb.org/download/{row.entry_pdbid.lower()}.pdb"
        request = requests.get(protein_link, '\n')
        with open(Path(output_folder) / f"{row.entry_pdbid.lower()}.pdb", 'w+') as f:
            f.write(request.text)
    protein = pd.parsePDB(str(pdb_file)).select(f"chain {row.chain} and protein")
    pdb_chain_filename = Path(output_folder) / f"{pdb_chain}.pdb"
    pd.writePDB(str(pdb_chain_filename), protein)
    return pdb_chain_filename


def smiles_to_3d_mol(smiles, output_sdf_file):
    """
    Convert SMILES representation to 3D molecule and write as SDF file
    """
    mol = Chem.MolFromSmiles(smiles)
    mol_H = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_H)

    with Chem.SDWriter(output_sdf_file) as w:
        w.write(mol_H)


def read_mol_from_sdf_file(sdf_file):
    """
    Obtain an RDKit molecule object from an SDF File
    """
    with Chem.SDMolSupplier(sdf_file) as r:
        return next(r)


def write_mol_to_sdf_file(mol, sdf_file):
    """
    Write an RDKit molecule object to an SDF File
    """
    with Chem.SDWriter(sdf_file) as w:
        w.write(mol)


def mol_to_smiles(sdf_file):
    """
    Convert an SDF file to a SMILES string
    """
    mol = read_mol_from_sdf_file(sdf_file)
    return Chem.MolToSmiles(mol)


def pdb_ligand_to_sdf(pdb_file, template_smiles, output_sdf_file):
    """
    Extract a ligand (with a given SMILES) from a PDB file as an SDF file
    """
    pdb = pd.parsePDB(str(pdb_file))
    ligand_residues = pdb.select("hetatm and not water")
    ligand_resnames = list(set(ligand_residues.getResnames()))
    template = AllChem.MolFromSmiles(template_smiles)
    rd_mol = None
    for name in ligand_resnames:
        # get pdb mol matching template mol
        ligand = pdb.select(f"resname {name}")
        output = StringIO()
        pd.writePDBStream(output, ligand)
        pdb_string = output.getvalue()
        rd_mol = AllChem.MolFromPDBBlock(pdb_string)
        if rd_mol.GetNumAtoms() == template.GetNumAtoms():
            break
    assert rd_mol is not None
    # remove stereochemistry
    for bond in rd_mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    Chem.SanitizeMol(rd_mol)
    RemoveStereochemistry(rd_mol)
    # assign bond orders from template mol
    new_mol = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
    AssignStereochemistryFrom3D(new_mol)
    # write to file
    with Chem.SDWriter(output_sdf_file) as w:
        w.write(new_mol)
        

def protein_ligand_csv(row, output_folder):
    """
    add path of pdb and sdf file to csv file
    """
    pdb_chain, ligand_resnum = row.PocketID.split(":")
    ligand_resnum = ligand_resnum.replace(".", "_")
    sdf_file = Path(output_folder) / f"{pdb_chain}__{ligand_resnum}.sdf"
    pdb_file = Path(output_folder) / f"{pdb_chain}.pdb"
    protein_ligand_file = Path(output_folder) / "protein_ligand.csv"

    if not protein_ligand_file.exists():
        with open(protein_ligand_file, 'w+') as c:
            writer = csv.writer(c)
            first_row = ['protein_path', 'ligand']
            writer.writerow(first_row)

    with open(protein_ligand_file, 'a') as f:
        writer = csv.writer(f)
        row_content = [pdb_file, sdf_file]
        writer.writerow(row_content)


def diffdock_csv(ref_sdf_file):
    """
    add path of pdb and sdf file to csv file
    """
    pdb_chain, ligand_resnum = ref_sdf_file.split("__")
    sdf_file = Path() / f"{ligand_resnum}"
    pdb_file = Path() / f"{pdb_chain}.pdb"
    protein_ligand_file = Path() / "protein_ligand.csv"

    if not protein_ligand_file.exists():
        with open(protein_ligand_file, 'w+') as c:
            writer = csv.writer(c)
            first_row = ['protein_path', 'ligand']
            writer.writerow(first_row)

    with open(protein_ligand_file, 'a') as f:
        writer = csv.writer(f)
        row_content = [pdb_file, sdf_file]
        writer.writerow(row_content)
    
