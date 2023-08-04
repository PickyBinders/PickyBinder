#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import glob
import csv
import re


def extract_vina_scores(vina_sdf_files):
    """
    Takes a list of vina sdf files and extracts the Autodock Vina scores.
    :return: list of lists for each vina sdf file containing 'Tool', 'Complex', 'Pocket', 'Rank',
             'Vina-free_energy', 'Vina-intermolecular_energy', 'Vina-internal_energy'
    """
    vina_scores = []
    for file in vina_sdf_files:
        with open(file) as f:
            pattern = "free_energy"
            for line in f:
                if re.search(pattern, line):
                    filename = file.split('/')[-1]
                    complex = filename.split('_pocket')[0]
                    pocket = filename.split('_')[-3]
                    rank = filename.split('_')[-1].split('.sdf')[0]
                    free_energy = line.split(':')[1].split(',')[0].strip()
                    intermolecular_energy = line.split(':')[2].split(',')[0].strip()
                    internal_energy = line.split(':')[3].split('}')[0].strip()

                    vina_scores.append(
                        ['vina', complex, pocket, rank, free_energy, intermolecular_energy, internal_energy]
                    )

    return vina_scores


def extract_smina_scores(smina_sdf_files):
    """
    Takes a list of smina sdf files and extracts the SMINA scores.
    :return: list of lists for each smina sdf file containing 'Tool', 'Complex', 'Pocket',
             'Rank', 'SMINA-minimizedAffinity'
    """
    smina_scores = []
    for file in smina_sdf_files:
        with open(file) as f:
            for line in (f.readlines()[-3:-2]):
                filename = file.split('/')[-1]
                complex = filename.split('_pocket')[0]
                pocket = filename.split('_')[-3]
                rank = filename.split('_')[-1].split('.sdf')[0]
                minimized_affinity = line.rstrip()

                smina_scores.append(
                    ['smina', complex, pocket, rank, minimized_affinity]
                )

    return smina_scores


def extract_gnina_scores(gnina_sdf_files):
    """
    Takes a list of gnina sdf files and extracts the GNINA scores.
    :return: list of lists for each gnina sdf file containing 'Tool', 'Complex', 'Pocket', 'Rank',
             'GNINA-minimizedAffinity', 'GNINA-CNNScore', 'GNINA-CNNAffinity'
    """
    gnina_scores = []
    for file in gnina_sdf_files:
        filename = file.split('/')[-1]
        complex = filename.split('_pocket')[0]
        pocket = filename.split('_')[-3]
        rank = filename.split('_')[-1].split('.sdf')[0]

        with open(file) as f:
            for line in (f.readlines()[-15:-14]):
                minimized_affinity = line.rstrip()

        with open(file) as f:
            for line in (f.readlines()[-12:-11]):
                cnn_score = line.rstrip()

        with open(file) as f:
            for line in (f.readlines()[-9:-8]):
                cnn_affinity = line.rstrip()

        gnina_scores.append(
            ['gnina', complex, pocket, rank, minimized_affinity, cnn_score, cnn_affinity]
        )

    return gnina_scores


def extract_tankbind_scores(tankbind_affinity_files):
    """
    Takes a list of tankbind affinity files and extracts the TANKBind affinity.
    :return: list of lists for each tankbind affinity file containing 'Tool', 'Complex', 'Pocket',
             'Rank', 'TANKBind-Affinity'
    """
    tb_scores = []
    for file in tankbind_affinity_files:
        with open(file) as f:
            csvreader = csv.reader(f)
            next(csvreader, None)   # skip the headers
            for index, row in enumerate(csvreader):
                tb_scores.append(
                    ['tankbind', row[1].split('_rdkit')[0], row[2], index + 1, row[4]]
                )

    return tb_scores


def extract_diffdock_scores(diffdock_sdf_names):
    """
    Takes a list of diffdock sdf file names and extracts the confidence score.
    :return: list of lists for each diffdock sdf file containing 'Tool', 'Complex', 'Rank', 'Confidence'
    """
    dd_scores = []
    for file in diffdock_sdf_names:
        complex_name = file.split('_rank')[0]
        rank = file.split('_rank')[1].split('_')[0]
        confidence = file.split('_confidence')[1].split('.sdf')[0]

        dd_scores.append(['diffdock', complex_name, rank, confidence])

    return dd_scores
