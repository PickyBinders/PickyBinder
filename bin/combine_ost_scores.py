#!/usr/bin/env python

import glob
import json
import os
import sys

import pandas as pd

tool_name = sys.argv[1]
complex_name = sys.argv[2]


def combine_scores(ost_json, tool):
    """
    """
    with open(ost_json, 'r') as f:
        data = json.load(f)

    name = data['model_ligands'][0]

    if tool == 'vina' or tool == 'smina' or tool == 'gnina':
        pocket = name.split('_')[-3]
        rank = name.split('_')[-1].split('.')[0]
        # confidence = '-'
    elif tool == 'tankbind':
        pocket = name.split('_')[-3] + name.split('_')[-2]
        rank = '-'
        # confidence = '-'
    elif tool == 'diffdock':
        pocket = '-'
        rank = name.split('_')[0].split('rank')[1]
        # confidence = name.split('-')[1].split('.sdf')[0]
    else:
        pocket = '-'
        rank = '-'

    if data['lddt_pli'] == {}:
        lddt_pli = '-'
        rmsd = '-'
        ref_ligand = ''
    else:
        lddt_pli = data['lddt_pli'][name]['lddt_pli']
        rmsd = data['lddt_pli'][name]['rmsd']
        ref_ligand = data['lddt_pli'][name]['reference_ligand']

    list_row = pd.DataFrame([[tool, complex_name, pocket, rank, lddt_pli, rmsd, ref_ligand]],
                            columns=['Tool', 'Complex', 'Pocket', 'Rank', 'lddt_pli', 'rmsd', 'Reference_Ligand'])
    # list_row = pd.DataFrame([[tool, complex_name, pocket, rank, confidence, lddt_pli, rmsd, ref_ligand]], columns=['Tool', 'Complex', 'Pocket', 'Rank', 'dd_confidence', 'lddt_pli', 'rmsd'])

    return list_row


files = [f for f in glob.glob("*.json")]
score_df = pd.DataFrame(columns=['Tool', 'Complex', 'Pocket', 'Rank', 'lddt_pli', 'rmsd', 'Reference_Ligand'])
# score_df = pd.DataFrame(columns=['Tool', 'Complex', 'Pocket', 'Rank', 'dd_confidence', 'lddt_pli', 'rmsd', 'Reference Ligand'])

for file in files:
    score_df = pd.concat([score_df, combine_scores(file, tool_name)])

score_df.Rank = pd.to_numeric(score_df.Rank, errors='coerce')
score_df.sort_values(by=['Pocket', 'Rank'], inplace=True)

if tool_name == 'vina' or tool_name == 'smina' or tool_name == 'gnina':
    pocket_nr = score_df.iloc[0]['Pocket']
    out_file = complex_name + '_' + pocket_nr + '_' + tool_name + '_score_summary.csv'
elif tool_name == 'tankbind' or tool_name == 'diffdock':
    out_file = complex_name + '_' + tool_name + '_score_summary.csv'

score_df.to_csv(out_file, index=False)
