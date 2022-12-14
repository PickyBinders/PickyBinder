#!/usr/bin/env python

from spyrmsd import io, rmsd
import os 
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# paths to the relevant data
ref_data = sys.argv[1]
diffdock_allData = sys.argv[2]
scoring_output = sys.argv[3]

# get all prediction folders
diffdock_predictions = [d for d in glob.glob(diffdock_allData + '/*/')]

"""rmsd prediction with sPyRMSD"""

for result in diffdock_predictions:
    # get all ranked ligand files
    files = [f for f in glob.glob(result + "*.sdf")]
    ligands = [os.path.split(path)[1] for path in files]
    ligands.sort()
    ligands.sort(key=len)

    # get the reference molecule
    normalized_path = os.path.normpath(result)
    path_components = normalized_path.split(os.sep)
    ref_mol = path_components[-1].split('-')[-1]

    # compare the structure of each predicted ligand against experimental data

    diffdock_dict = {}

    for ligand in ligands:
        diffdock_path = os.path.join(result, ligand)

        ref = io.loadmol(ref_data + '/' + ref_mol)
        ref.strip()

        diffdock = io.loadmol(diffdock_path)
        diffdock.strip()

        coords_ref = ref.coordinates
        anum_ref = ref.atomicnums
        adj_ref = ref.adjacency_matrix

        coords_diffdock = diffdock.coordinates
        anum_diffdock = diffdock.atomicnums
        adj_diffdock = diffdock.adjacency_matrix

        RMSD_diffdock = rmsd.symmrmsd(
            coords_ref,
            coords_diffdock,
            anum_ref,
            anum_diffdock,
            adj_ref,
            adj_diffdock,
        )

        # add confidence and RMSD value to the dict
        confidence = ligand.split('.sdf')[0].split('confidence')[-1]
        scoring = [confidence, RMSD_diffdock]
        diffdock_dict[ligand.split('_')[0]] = scoring

    # create table for the prediction
    diffdock_table = pd.DataFrame.from_dict(diffdock_dict, orient='index', columns=['Confidence', 'DiffDock RMSD'])
    diffdock_table.reset_index(inplace=True)
    diffdock_table = diffdock_table.rename(columns={'index': 'Diffdock rank'})
    # remove duplicated rank1
    diffdock_table = diffdock_table.iloc[1:, :]

    # save table as csv file
    diffdock_table.to_csv(scoring_output + '/' + ref_mol.split('.')[0] + '.csv', index=False)

    """add sucess rate"""
    # add columns for success interpretation
    diffdock_table['DiffDock docking successful'] = diffdock_table['DiffDock RMSD'] < 2

    # get metrics of successfull docking
    diffdock_success = diffdock_table['DiffDock docking successful'][
        diffdock_table['DiffDock docking successful'] == True].count()
    diffdock_no_success = diffdock_table['DiffDock docking successful'][
        diffdock_table['DiffDock docking successful'] == False].count()

    diffdock_success_perc = np.round(diffdock_success / (diffdock_success + diffdock_no_success) * 100, 1)
    diffdock_no_success_perc = np.round(diffdock_no_success / (diffdock_success + diffdock_no_success) * 100, 1)

    docking_success = pd.DataFrame(
        {'Tool': ['DiffDock'], 'Success': [diffdock_success], 'No Success': [diffdock_no_success],
         'Success Percent': [diffdock_success_perc], 'No Success Percent': [diffdock_no_success_perc]})

    """plots for RMSD overview and success rates of the tools"""
    # RMSDs
    plt.plot(diffdock_table['Diffdock rank'], diffdock_table['DiffDock RMSD'], 'bo', label='DiffDock')
    plt.ylabel = 'RMDS[Å]'
    plt.axis([-1, 41, 0, 30])
    plt.grid(True)
    plt.axhline(y=2, color='red', linestyle='--')
    plt.tick_params(axis='x', labelrotation=90)
    plt.legend(loc='upper center', ncol=2, shadow=True)
    plt.title('RMSD of DiffDock against experimental data for ' + ref_mol.split('.')[0])
    plt.savefig(scoring_output + '/' + ref_mol.split('.')[0] + '_RMSD.png')
    plt.close()

    # Performance overview
    plt.barh(docking_success['Tool'], docking_success['Success Percent'], color='green')
    plt.barh(docking_success['Tool'], docking_success['No Success Percent'], left=docking_success['Success Percent'],
             color='grey')
    plt.legend(['RMDS[Å] < 2', 'RMDS[Å] > 2'], bbox_to_anchor=(0.75, 0.5), ncol=2)
    plt.title('Perfomrance of Diffdock for ' + ref_mol.split('.')[0])
    plt.savefig(scoring_output + '/' + ref_mol.split('.')[0] + '_performance.png')
    plt.close()

"""create summary table for all molecules"""

csv_tables = sys.argv[3]

summary_table = pd.DataFrame(
    columns=['Molecule', 'DiffDock rank1 RMSD', 'Diffdock best RMSD', 'Diffdock best RMSD rank'])

# get the paths of the csv files
files = [f for f in glob.glob(csv_tables + '/' + "*.csv")]
if csv_tables + '/' +'summary_table.csv' in files:
    files.remove(csv_tables + '/' + 'summary_table.csv')

for table in files:
    df = pd.read_csv(table)
    df1 = df.drop(['Confidence'], axis=1)

    molecule = table.split('/')[-1].split('.')[0]

    df_sorted = df1.sort_values(by='DiffDock RMSD', ignore_index=True)
    best = df_sorted.loc[0, 'DiffDock RMSD']
    best_rank = df_sorted.loc[0, 'Diffdock rank']

    df2 = df1.loc[df['Diffdock rank'] == 'rank1']
    summary = df2.copy()
    summary.rename(columns={'Diffdock rank': 'Molecule', 'DiffDock RMSD': 'DiffDock rank1 RMSD'}, inplace=True)
    summary.loc[0, 'Molecule'] = molecule

    summary.loc[0, 'Diffdock best RMSD'] = best
    summary.loc[0, 'Diffdock best RMSD rank'] = best_rank

    summary_table = pd.concat([summary_table, summary])

# save summary table as csv
summary_table.to_csv(scoring_output + '/' + 'summary_table.csv', index=False)
