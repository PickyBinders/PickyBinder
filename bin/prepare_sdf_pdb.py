from ligand_utility import get_ligand_sdf, get_protein_pdb, mol_to_smiles, protein_ligand_csv
import pandas as pd
import sys

dataset = sys.argv[1]
output_folder = sys.argv[2]

df = pd.read_csv(dataset)


for row in range(len(df.index)):

    get_ligand_sdf(df.loc[df.index[row]], output_folder)

    get_protein_pdb(df.loc[df.index[row]], output_folder)

    protein_ligand_csv(df.loc[df.index[row]], output_folder)
