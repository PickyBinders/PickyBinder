#!/usr/bin/env python

from ligand_utility import diffdock_csv
from pathlib import Path
import sys

naming = sys.argv[1]
receptor_hs = sys.argv[2]
ref_sdf_files = sys.argv[3:]

for ref in ref_sdf_files:
    diffdock_csv(ref, naming, receptor_hs)
