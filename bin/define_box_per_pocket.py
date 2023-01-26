#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import sys

complex = sys.argv[1]
p2rank_predictions = sys.argv[2]
box_size = sys.argv[3]

with open(box_size, 'r') as infile:
  size = infile.read().strip()
  
  size_x = "size_x = " + size
  size_y = "size_y = " + size
  size_z = "size_z = " + size


df = pd.read_csv(p2rank_predictions)
df.columns = df.columns.str.strip()
df["name"] = df["name"].str.strip()

for row in range(len(df.index)):
    pocket = df.loc[row,"name"]
    pocket_file = Path() / f"{complex}_{pocket}.box"
    print(pocket_file)
    
    center_x = 'center_x = ' + str(df.loc[row,"center_x"])
    center_y = 'center_y = ' + str(df.loc[row,"center_y"])
    center_z = 'center_z = ' + str(df.loc[row,"center_z"])

    box = [center_x, center_y, center_z, size_x, size_y, size_z]
    
    with open(pocket_file, 'w+') as f:
        f.write('\n'.join(box))
        f.write("\n")
