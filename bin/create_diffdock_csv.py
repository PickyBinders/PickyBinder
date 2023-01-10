#!/usr/bin/env python

from ligand_utility import diffdock_csv
from pathlib import Path
import sys


ref_sdf_files = sys.argv[1:]

for ref in ref_sdf_files:
  diffdock_csv(ref)
