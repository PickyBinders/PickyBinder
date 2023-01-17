#!/usr/bin/env python

from pathlib import Path
import prody as pd
import numpy as np
import sys

receptor_chain = sys.argv[1]
 
protein, prot_chain = receptor_chain.split('_')
protein = protein.lower()

atoms = pd.parsePDB(protein, chain=prot_chain)
center = pd.calcCenter(atoms).round(3)

center_x = 'center_x = ' + str(center[0])
center_y = 'center_y = ' + str(center[1])
center_z = 'center_z = ' + str(center[2])

lines = [center_x, center_y, center_z, "size_x = 40.0", "size_y = 40.0", "size_z = 40.0", ""]

box_file = Path() / f"{receptor_chain}_box.txt"

with open(box_file, 'w+') as f:
  f.write('\n'.join(lines))
