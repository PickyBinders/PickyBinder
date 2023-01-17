#!/usr/bin/env python

from ligand_utility import get_ligand_sdf, get_protein_pdb, mol_to_smiles, protein_ligand_csv
import pandas as pd
import sys

dataset = sys.argv[1]
output_folder = sys.argv[2]

df = pd.read_csv(dataset)

for row in range(len(df.index)):
    
    try:
      get_ligand_sdf(df.loc[row], output_folder)
    except (AttributeError, TypeError) as err:
      print("get_ligand_sdf: ", df.loc[row].iat[0])
      print(err)
    
    try:
      get_protein_pdb(df.loc[row], output_folder)
    except (AttributeError, TypeError) as err:
      print("get_protein_pdb: ", df.loc[row].iat[0])
      print(err)
    
    try:
      protein_ligand_csv(df.loc[row], output_folder)
    except (AttributeError, TypeError) as err:
      print("protein_ligand_csv: ", df.loc[row].iat[0])
      print(err)
