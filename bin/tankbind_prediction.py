#!/usr/bin/env python

"""
This script is adapted from TankBind
https://github.com/luwei0917/TankBind/blob/main/examples/prediction_example_using_PDB_6hd6.ipynb
as of 06.02.2023 
"""

import sys
import os
import numpy as np
import pandas as pd
import torch
import logging

import rdkit.Chem as Chem
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser
from torch_geometric.loader import DataLoader
from tqdm import tqdm

tankbind_scripts = sys.argv[1]
sys.path.insert(0, tankbind_scripts)

from feature_utils import get_protein_feature, get_clean_res_list, extract_torchdrug_feature_from_mol
from data import TankBind_prediction
from model import get_model
from generation_utils import get_LAS_distance_constraint_mask, get_info_pred_distance, write_with_new_coords

# renaming input
complex_name = sys.argv[2]
protein_file = sys.argv[3]
ligand_file = sys.argv[4]
p2rank_predictions = sys.argv[5]
naming = sys.argv[6]

if naming == 'default':
    pdb, ligand = complex_name.split('__')
else:
    pdb = complex_name.split('_')[0]
    ligand = complex_name

# get protein feature
parser = PDBParser(QUIET=True)
s = parser.get_structure("x", protein_file)
res_list = get_clean_res_list(s.get_residues(), verbose=False, ensure_ca_exist=True)

protein_dict = {}
protein_dict[pdb] = get_protein_feature(res_list)

# get compound feature
compound_dict = {}
mol = Chem.MolFromMolFile(ligand_file)
mol = Chem.RemoveHs(mol)
compound_dict[pdb + f"_{ligand}" + "_rdkit"] = extract_torchdrug_feature_from_mol(mol, has_LAS_mask=True)

# p2rank info
info = []
for compound_name in list(compound_dict.keys()):
    # use protein center as the block center.
    com = ",".join([str(a.round(3)) for a in protein_dict[pdb][0].mean(axis=0).numpy()])
    info.append([pdb, compound_name, "protein_center", com])

    pocket = pd.read_csv(p2rank_predictions)
    pocket.columns = pocket.columns.str.strip()
    pocket_coms = pocket[['center_x', 'center_y', 'center_z']].values
    for ith_pocket, com in enumerate(pocket_coms):
        com = ",".join([str(a.round(3)) for a in com])
        info.append([pdb, compound_name, f"pocket_{ith_pocket + 1}", com])
info = pd.DataFrame(info, columns=['protein_name', 'compound_name', 'pocket_name', 'pocket_com'])

# construct dataset

# check how many threads possible/useful
torch.set_num_threads(1)

dataset_path = f"{complex_name}_tankbindDataset/"
os.system(f"mkdir -p {dataset_path}")
dataset = TankBind_prediction(dataset_path, data=info, protein_dict=protein_dict, compound_dict=compound_dict)

batch_size = 5
# device = 'cuda' if torch.cuda.is_available() else 'cpu'
device = 'cpu'
logging.basicConfig(level=logging.INFO)
model = get_model(0, logging, device)

# self-dock model
modelFile = f"{tankbind_scripts}/saved_models/self_dock.pt"

# re-dock model
# modelFile = f"{tankbind_scripts}/saved_models/re_dock.pt"

model.load_state_dict(torch.load(modelFile, map_location=device))
_ = model.eval()

# if problem occur with DataLoader change num_workers to 0
data_loader = DataLoader(dataset, batch_size=batch_size, follow_batch=['x', 'y', 'compound_pair'], shuffle=False,
                         num_workers=0)
affinity_pred_list = []
y_pred_list = []
for data in tqdm(data_loader):
    data = data.to(device)
    y_pred, affinity_pred = model(data)
    affinity_pred_list.append(affinity_pred.detach().cpu())
    for i in range(data.y_batch.max() + 1):
        y_pred_list.append((y_pred[data['y_batch'] == i]).detach().cpu())

affinity_pred_list = torch.cat(affinity_pred_list)

info = dataset.data
info['affinity'] = affinity_pred_list
info_sorted = info.sort_values('affinity', ascending=False)

info_sorted.to_csv(f"{complex_name}_tankbind.csv")


# from predicted interaction distance map to sdf
device = 'cpu'
for idx, (dataframe_idx, line) in enumerate(info_sorted.iterrows()):
    pocket_name = line['pocket_name']
    compound_name = line['compound_name']
    ligandName = compound_name.split("_")[1]
    coords = dataset[dataframe_idx].coords.to(device)
    protein_nodes_xyz = dataset[dataframe_idx].node_xyz.to(device)
    n_compound = coords.shape[0]
    n_protein = protein_nodes_xyz.shape[0]
    y_pred = y_pred_list[dataframe_idx].reshape(n_protein, n_compound).to(device)
    y = dataset[dataframe_idx].dis_map.reshape(n_protein, n_compound).to(device)
    compound_pair_dis_constraint = torch.cdist(coords, coords)
    mol = Chem.MolFromMolFile(ligand_file)
    mol = Chem.RemoveHs(mol)
    LAS_distance_constraint_mask = get_LAS_distance_constraint_mask(mol).bool()
    pred_info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz, compound_pair_dis_constraint,
                                       LAS_distance_constraint_mask=LAS_distance_constraint_mask,
                                       n_repeat=1, show_progress=False)

    result_folder = 'tankbind_predictions/'
    os.system(f'mkdir -p {result_folder}')
    toFile = f'{result_folder}/{complex_name}_{pocket_name}_tankbind.sdf'
    new_coords = pred_info.sort_values("loss")['coords'].iloc[0].astype(np.double)
    write_with_new_coords(mol, new_coords, toFile)
