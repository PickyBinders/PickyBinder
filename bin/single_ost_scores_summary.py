#!/usr/bin/env python

import glob
import json
import os
import sys

import pandas as pd


def combine_scores(ost_json):
    """
    """
    with open(ost_json, 'r') as f:
        data = json.load(f)

    try:
        name = data['model_ligands'][0]
        complex_name_parts = ost_json.split('.json')[0].split('_')
        complex_name = '_'.join(complex_name_parts[:-1])
        rank = complex_name_parts[-1]

        if data['lddt_pli'] == {}:
            lddt_pli = '-'
            lddt_lp = '-'
            rmsd = '-'
            ref_ligand = ''
        else:
            lddt_pli = data['lddt_pli'][name]['lddt_pli']
            lddt_lp = data['lddt_pli'][name]['lddt_lp']
            rmsd = data['lddt_pli'][name]['rmsd']
            ref_ligand = data['lddt_pli'][name]['reference_ligand']

        list_row = pd.DataFrame([[complex_name, rank, lddt_pli, lddt_lp, rmsd, ref_ligand]],
                                columns=['Complex', 'Rank', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand'])

        return list_row
    except:
        pass


files = [f for f in glob.glob("*.json")]

score_df = pd.DataFrame(columns=['Complex', 'Rank', 'lDDT-PLI', 'lDDT-LP', 'BiSyRMSD', 'Reference_Ligand'])

for file in files:
    score_df = pd.concat([score_df, combine_scores(file)])

score_df.Rank = pd.to_numeric(score_df.Rank, errors='coerce')
score_df.sort_values(by=['Complex', 'Rank'], inplace=True)

out_file = 'score_summary.csv'

score_df.to_csv(out_file, index=False)
