#!/usr/bin/env python

import csv
import glob
import os
from pathlib import Path

import pandas as pd
from get_tool_scores import (
    extract_tankbind_scores,
    extract_diffdock_scores,
    extract_vina_scores,
    extract_smina_scores,
    extract_gnina_scores,
)


# Set input and output file paths
ligand_score_summary_file = 'ligand_score_summary.csv'
out_file = 'all_scores_summary.csv'

# Load the ligand score summary file
all_scores = pd.read_csv(ligand_score_summary_file)

# Convert columns to appropriate types
all_scores['Rank'] = all_scores['Rank'].fillna(0).astype(int).astype(str)
all_scores['lddt_pli'] = all_scores['lddt_pli'].astype(str)
all_scores['rmsd'] = all_scores['rmsd'].astype(str)
all_scores['center_x'] = all_scores['center_x'].astype(str)
all_scores['center_y'] = all_scores['center_y'].astype(str)
all_scores['center_z'] = all_scores['center_z'].astype(str)

# Tankbind
tankbind_files = [f for f in glob.glob("*_tankbind.csv")]

if tankbind_files:
    tankbind_scores_list = extract_tankbind_scores(tankbind_files)
    tankbind_scores = pd.DataFrame(
        tankbind_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'TANKBind-Affinity']
    )
    all_scores = tankbind_scores.set_index(['Tool', 'Complex', 'Pocket']).combine_first(
        all_scores.set_index(['Tool', 'Complex', 'Pocket'])).reset_index()
    all_scores = all_scores.reindex(
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'lddt_pli', 'rmsd', 'Reference_Ligand', 'center_x', 'center_y',
                 'center_z', 'TANKBind-Affinity']
    )

# Diffdock
diffdock_file = 'dd_file_names.txt'

if os.path.exists(diffdock_file):
    with open(diffdock_file) as f:
        diffdock_files = [line.rstrip() for line in f]
    diffdock_scores_list = extract_diffdock_scores(diffdock_files)
    diffdock_scores = pd.DataFrame(
        diffdock_scores_list,
        columns=['Tool', 'Complex', 'Rank', 'DiffDock-Confidence']
    )
    all_scores = pd.merge(all_scores, diffdock_scores, how="outer", on=["Tool", "Complex", "Rank"])

# Vina
vina_sdfs = [f for f in glob.glob("*_vina_*.sdf")]
if vina_sdfs:
    vina_scores_list = extract_vina_scores(vina_sdfs)
    vina_scores = pd.DataFrame(
        vina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'Vina-free_energy', 'Vina-intermolecular_energy',
                 'Vina-internal_energy']
    )
    all_scores = pd.merge(all_scores, vina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])

# SMINA
smina_sdfs = [f for f in glob.glob("*_smina_*.sdf")]
if smina_sdfs:
    smina_scores_list = extract_smina_scores(smina_sdfs)
    smina_scores = pd.DataFrame(
        smina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'SMINA-minimizedAffinity']
    )
    all_scores = pd.merge(all_scores, smina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])

# GNINA
gnina_sdfs = [f for f in glob.glob("*_gnina_*.sdf")]
if gnina_sdfs:
    gnina_scores_list = extract_gnina_scores(gnina_sdfs)
    gnina_scores = pd.DataFrame(
        gnina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'GNINA-minimizedAffinity', 'GNINA-CNNScore', 'GNINA-CNNAffinity']
    )
    all_scores = pd.merge(all_scores, gnina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])

# Write the combined scores to the output file

all_scores.to_csv(out_file, index=False)
