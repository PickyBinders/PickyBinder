#!/usr/bin/env python

import csv
import glob
import os
import sys
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
launchdir = sys.argv[1]
ligand_score_summary_file = 'ligand_score_summary.csv'

# Load the ligand score summary file
ost_scores = pd.read_csv(ligand_score_summary_file)

# Convert columns to appropriate types
ost_scores['Rank'] = ost_scores['Rank'].fillna(0).astype(int).astype(str)
ost_scores['lDDT-PLI'] = ost_scores['lDDT-PLI'].astype(str)
ost_scores['lDDT-LP'] = ost_scores['lDDT-LP'].astype(str)
ost_scores['BiSyRMSD'] = ost_scores['BiSyRMSD'].astype(str)
ost_scores['Box_Center_x'] = ost_scores['Box_Center_x'].astype(str)
ost_scores['Box_Center_y'] = ost_scores['Box_Center_y'].astype(str)
ost_scores['Box_Center_z'] = ost_scores['Box_Center_z'].astype(str)

# Tankbind
tankbind_files = [f for f in glob.glob(launchdir + "/predictions/tankbind/*/*_tankbind.csv")]

if tankbind_files:
    tankbind_scores_list = extract_tankbind_scores(tankbind_files)
    tankbind_scores = pd.DataFrame(
        tankbind_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'TANKBind-Affinity']
    )
    tb_summary = ost_scores[ost_scores.Tool == 'tankbind'].copy()
    tb_summary = tankbind_scores.set_index(['Tool', 'Complex', 'Pocket']).combine_first(
        tb_summary.set_index(['Tool', 'Complex', 'Pocket'])).reset_index()
    tb_summary = tb_summary.reindex(
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand', 'Box_Center_x',
                 'Box_Center_y', 'Box_Center_z', 'TANKBind-Affinity']
    )
    tb_summary = tb_summary[['Tool', 'Complex', 'Pocket', 'Box_Center_x', 'Box_Center_y', 'Box_Center_z', 'Rank',
                             'TANKBind-Affinity', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand']]
    tb_summary.to_csv('tankbind_summary.csv', index=False)

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
    dd_summary = ost_scores[ost_scores.Tool == 'diffdock'].copy()
    dd_summary = pd.merge(dd_summary, diffdock_scores, how="outer", on=["Tool", "Complex", "Rank"])
    dd_summary = dd_summary[['Tool', 'Complex', 'Rank', 'DiffDock-Confidence', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD',
                             'Reference_Ligand']]
    dd_summary.to_csv('diffdock_summary.csv', index=False)

# Vina
vina_sdfs = [f for f in glob.glob(launchdir + "/predictions/vina/vina_predictions/*/*/*_vina_*.sdf")]
if vina_sdfs:
    vina_scores_list = extract_vina_scores(vina_sdfs)
    vina_scores = pd.DataFrame(
        vina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'Vina-free_energy', 'Vina-intermolecular_energy',
                 'Vina-internal_energy']
    )
    vina_summary = ost_scores[ost_scores.Tool == 'vina'].copy()
    vina_summary = pd.merge(vina_summary, vina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])
    vina_summary = vina_summary[['Tool', 'Complex', 'Pocket', 'Box_Center_x', 'Box_Center_y', 'Box_Center_z', 'Rank',
                                 'Vina-free_energy', 'Vina-intermolecular_energy', 'Vina-internal_energy', 'lDDT-PLI',
                                 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand']]
    vina_summary.to_csv('vina_summary.csv', index=False)

# SMINA
smina_sdfs = [f for f in glob.glob(launchdir + "/predictions/smina/*/*/*_smina_*.sdf")]
if smina_sdfs:
    smina_scores_list = extract_smina_scores(smina_sdfs)
    smina_scores = pd.DataFrame(
        smina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'SMINA-minimizedAffinity']
    )
    smina_summary = ost_scores[ost_scores.Tool == 'smina'].copy()
    smina_summary = pd.merge(smina_summary, smina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])
    smina_summary = smina_summary[['Tool', 'Complex', 'Pocket', 'Box_Center_x', 'Box_Center_y', 'Box_Center_z', 'Rank',
                                   'SMINA-minimizedAffinity', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand']]
    smina_summary.to_csv('smina_summary.csv', index=False)

# GNINA
gnina_sdfs = [f for f in glob.glob(launchdir + "/predictions/gnina/*/*/*_gnina_*.sdf")]
if gnina_sdfs:
    gnina_scores_list = extract_gnina_scores(gnina_sdfs)
    gnina_scores = pd.DataFrame(
        gnina_scores_list,
        columns=['Tool', 'Complex', 'Pocket', 'Rank', 'GNINA-minimizedAffinity', 'GNINA-CNNScore', 'GNINA-CNNAffinity']
    )
    gnina_summary = ost_scores[ost_scores.Tool == 'gnina'].copy()
    gnina_summary = pd.merge(gnina_summary, gnina_scores, how="outer", on=["Tool", "Complex", "Pocket", "Rank"])
    gnina_summary = gnina_summary[['Tool', 'Complex', 'Pocket', 'Box_Center_x', 'Box_Center_y', 'Box_Center_z', 'Rank',
                                   'GNINA-minimizedAffinity', 'GNINA-CNNScore', 'GNINA-CNNAffinity', 'lDDT-PLI',
                                   'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand']]
    gnina_summary.to_csv('gnina_summary.csv', index=False)
