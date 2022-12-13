from ligand_utility import get_ligand_sdf, get_protein_pdb, mol_to_smiles, protein_ligand_csv
import pandas as pd

dataset = '/Users/mleemann/Documents/test_data/df_entry_ligand_pocket_edia_10.csv'

df = pd.read_csv(dataset)
output_folder = '/Users/mleemann/Documents/test_data/test_output'

for row in range(len(df.index)):

    get_ligand_sdf(df.loc[df.index[row]], output_folder)

    get_protein_pdb(df.loc[df.index[row]], output_folder)

    protein_ligand_csv(df.loc[df.index[row]], output_folder)